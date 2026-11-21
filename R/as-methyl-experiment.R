#' Align reference methylation calls across samples
#'
#' Project regional calls onto a shared set of reference sites. DRAGEN inputs
#' contribute pooled reference-allele measurements. Rastair inputs contribute
#' reference CpG context measurements. All inputs must use the same format
#' and reference genome.
#'
#' @param calls A `GRanges` from [readMethylVcf()] containing one or more VCF
#'   samples, or a named list of these objects with one sample per element.
#'   List names become column names. Use the reader's `samples` argument to
#'   select a single sample from a multisample file for a list element.
#' @param filters Exact VCF FILTER values to retain. By default, retain
#'   `"PASS"` and `"."` (filters not applied). Use `NULL` to include all
#'   records.
#'
#' @return A [SummarizedExperiment::RangedSummarizedExperiment] with one row
#'   per chromosome, 1-based anchor, strand, context, and modification across
#'   the union of selected sites. Each selected sample has a column, including
#'   samples without calls. `colData()` stores the original VCF sample ID and
#'   an `import_index` into `S4Vectors::metadata(result)$imports`, which
#'   retains the original imports' file metadata and diagnostics. The
#'   experiment metadata also record the chosen `filters`.
#'
#' @details The `fraction` assay contains fractions backed by positive
#'   informative depth. The `reported_fraction` assay preserves the raw VCF
#'   values. `informative_depth`, `modified_depth` and `coverage_status`
#'   preserve the original values independently. All assays use `NA` where
#'   there is no call for a site and sample. A reported zero-depth call keeps
#'   its depth of zero and its coverage status.
#'
#'   Strand-specific calls and collapsed dyads occupy separate rows. This
#'   projection excludes alternate, de novo, allele-specific and haplotype
#'   measurements. Use the original ranges to analyze those measurements.
#'   Multiple calls that pass the `filters` selection at one site in one
#'   sample signal an error with class `methylVCF_duplicate_site`.
#'
#'   Conflicting assemblies or sequence lengths signal an error. Different
#'   reference file locations are allowed when the inputs declare a common
#'   assembly. Other conflicting `##reference` declarations signal an error.
#'   Confirm a common assembly when header information is incomplete. The
#'   assays use dense matrices intended for regional data.
#'
#' @examples
#' file <- system.file(
#'     "extdata", "dragen-multisample.g.vcf", package = "methylVCF"
#' )
#' region <- GenomicRanges::GRanges("chr21", IRanges::IRanges(13204474L))
#' calls <- readMethylVcf(file, region)
#' experiment <- asMethylExperiment(calls)
#' SummarizedExperiment::assay(experiment, "fraction")
#' SummarizedExperiment::assay(experiment, "reported_fraction")
#'
#' @testexamples
#' testthat::expect_s4_class(experiment, "RangedSummarizedExperiment")
#' testthat::expect_identical(dim(experiment), c(1L, 2L))
#' testthat::expect_identical(colnames(experiment), c("sample_a", "sample_b"))
#' testthat::expect_equal(
#'     unname(SummarizedExperiment::assay(experiment, "fraction")[1L, ]),
#'     c(0.2381, NA_real_)
#' )
#' testthat::expect_equal(
#'     unname(SummarizedExperiment::assay(experiment, "reported_fraction")[1L, ]),
#'     c(0.2381, 0)
#' )
#'
#' @seealso [readMethylVcf()] for importing the individual measurements and
#'   `vignette("methylVCF")` for the complete workflow.
#' @export
#' @md
asMethylExperiment <- function(calls, filters = c("PASS", ".")) {
    validate_methyl_filters(filters)
    single <- inherits(calls, "GRanges")
    inputs <- if (single) list(calls) else calls
    valid_inputs <- is.list(inputs) && length(inputs) > 0L && all(vapply(inputs, inherits, logical(1L), "GRanges"))
    if (!valid_inputs) {
        cli::cli_abort(
            "{.arg calls} must be a GRanges or a named list of GRanges.",
            call = NULL
        )
    }
    if (!single && !is_valid_sample_ids(names(inputs))) {
        cli::cli_abort(
            "The list in {.arg calls} needs unique, nonempty sample names.",
            call = NULL
        )
    }

    imports <- lapply(inputs, S4Vectors::metadata)
    formats <- lapply(imports, function(info) info[["format"]])
    format <- unique(unlist(formats, use.names = FALSE))
    same_format <- all(lengths(formats) == 1L) && length(format) == 1L && format %in% c("dragen", "rastair")
    if (!same_format) {
        cli::cli_abort(
            "Use readMethylVcf() calls from one format, DRAGEN or Rastair.",
            call = NULL
        )
    }
    sequence_info <- merge_experiment_seqinfo(inputs, imports)
    sample_data <- experiment_sample_data(inputs, single)
    samples <- rownames(sample_data)
    sample_override <- if (single) list(NULL) else as.list(samples)
    scope <- if (format == "dragen") "pooled" else "context"
    rows <- Map(
        reference_measurements, 
        inputs, 
        sample_override,
        MoreArgs = list(scope = scope, filters = filters)
    )
    measurements <- do.call(rbind, unname(rows))
    sites <- experiment_sites(measurements, sequence_info)
    SummarizedExperiment::SummarizedExperiment(
        assays = experiment_assays(measurements, sites, samples),
        rowRanges = sites,
        colData = sample_data,
        metadata = list(
            format = format, 
            filters = filters, 
            imports = unname(imports)
        )
    )
}

#' Map experiment columns back to the sample selected in each VCF import
#'
#' @param inputs A list of regional calls with VCF sample metadata.
#' @param single Whether the user supplied one multisample query directly.
#' @return A DataFrame in assay column order with VCF IDs and import indices.
#' @noRd
experiment_sample_data <- function(inputs, single) {
    ids <- lapply(inputs, experiment_sample_names)
    if (single) {
        samples <- ids[[1L]]
        original_ids <- samples
        import_index <- rep(1L, length(samples))
    } else {
        if (any(lengths(ids) > 1L)) {
            cli::cli_abort(
                "Select one VCF sample per list element in {.arg calls}.",
                call = NULL
            )
        }
        samples <- names(inputs)
        original_ids <- vapply(ids, `[`, character(1L), 1L)
        import_index <- seq_along(inputs)
    }
    S4Vectors::DataFrame(
        vcf_sample = original_ids, 
        import_index = import_index,
        row.names = samples
    )
}

#' Combine the reference information from a set of regional queries
#'
#' @param inputs List of input `GRanges`.
#' @param imports File metadata for each input.
#' @return A merged `Seqinfo`, or an error for incompatible declarations.
#' @noRd
merge_experiment_seqinfo <- function(inputs, imports) {
    sequence_info <- lapply(inputs, GenomicRanges::seqinfo)
    assemblies <- lapply(sequence_info, function(info) {
        labels <- Seqinfo::genome(info)
        unique(labels[!is.na(labels) & nzchar(labels)])
    })
    known <- lengths(assemblies) > 0L
    if (sum(known) > 1L && !length(Reduce(intersect, assemblies[known]))) {
        labels <- unique(unlist(assemblies, use.names = FALSE))
        cli::cli_abort(
            "The inputs declare incompatible genome assemblies: {.val {labels}}.",
            call = NULL
        )
    }

    references <- vapply(imports, function(info) {
        info[["reference"]] %||% NA_character_
    }, character(1L))
    declared <- !is.na(references) & nzchar(references)
    distinct_references <- unique(references[declared])
    locations <- grepl("[/\\\\]", references[declared])
    compatible_locations <- all(locations) && all(known[declared])
    if (length(distinct_references) > 1L && !compatible_locations) {
        cli::cli_abort(
            "The inputs declare different VCF references: {.val {distinct_references}}.",
            call = NULL
        )
    }
    tryCatch(
        Reduce(merge, sequence_info),
        error = function(error) {
            cli::cli_abort(
                "The inputs have incompatible sequence information.",
                parent = error, call = NULL
            )
        }
    )
}

#' Recover the selected VCF samples including those without measurements
#'
#' @param input A single set of calls from the regional reader.
#' @return VCF sample IDs in query order, with a fallback to IDs present in
#'   the rows for results created before sample metadata was recorded.
#' @noRd
experiment_sample_names <- function(input) {
    if (!"sample" %in% names(S4Vectors::mcols(input))) {
        cli::cli_abort(
            "Each input must contain the sample column from readMethylVcf().",
            call = NULL
        )
    }
    observed <- unique(input$sample)
    samples <- S4Vectors::metadata(input)[["samples"]] %||% observed
    if (!is_valid_sample_ids(samples) || !all(observed %in% samples)) {
        cli::cli_abort(
            "The input contains invalid or inconsistent VCF sample IDs.",
            call = NULL
        )
    }
    samples
}

#' Select the comparable reference measurements for one import
#'
#' @param input A single set of calls from the regional reader.
#' @param sample Sample ID for a list element, or `NULL` to retain the VCF IDs.
#' @param scope Eligible evidence scope for this format.
#' @param filters Eligible exact VCF FILTER values, or `NULL` for all records.
#' @return A data frame with site coordinates, sample IDs and assay values.
#' @noRd
reference_measurements <- function(input, sample, scope, filters) {
    columns <- c(
        "sample", "scope", "motif_origin", "context", "modification",
        "fraction", "informative_depth", "modified_depth", "coverage_status"
    )
    fields <- S4Vectors::mcols(input)
    has_columns <- all(c(columns, "filter") %in% names(fields))
    base_anchors <- GenomicRanges::width(input) == 1L & GenomicRanges::start(input) >= 1L
    if (!has_columns || !all(base_anchors)) {
        cli::cli_abort(
            "Each input must contain one-base calls from readMethylVcf().",
            call = NULL
        )
    }
    eligible <- fields$scope == scope & fields$motif_origin == "reference"
    if (!is.null(filters)) { eligible <- eligible & fields$filter %in% filters }
    selected <- which(eligible)
    input <- input[selected]
    values <- as.data.frame(fields[selected, columns, drop = FALSE])
    if (!is.null(sample)) { values$sample <- rep(sample, nrow(values)) }
    data.frame(
        seqname = as.character(GenomicRanges::seqnames(input)),
        position = GenomicRanges::start(input),
        strand = as.character(GenomicRanges::strand(input)),
        values[setdiff(columns, c("scope", "motif_origin"))]
    )
}

#' Build one named genomic range per reference site
#'
#' @param measurements Combined reference calls.
#' @param sequence_info Merged sequence information for the imports.
#' @return A `GRanges` ordered by contig, anchor, strand and context, with
#'   row names unique for each site and modification.
#' @noRd
experiment_sites <- function(measurements, sequence_info) {
    keys <- c("seqname", "position", "strand", "context", "modification")
    sites <- unique(measurements[keys])
    contig <- match(sites$seqname, Seqinfo::seqlevels(sequence_info))
    site_order <- order(
        contig, 
        sites$position, 
        sites$strand, 
        sites$context, 
        sites$modification
    )
    sites <- sites[site_order, ]
    ranges <- GenomicRanges::GRanges(
        sites$seqname,
        IRanges::IRanges(sites$position, width = 1L),
        strand = sites$strand,
        seqinfo = sequence_info,
        context = sites$context,
        modification = sites$modification
    )
    names(ranges) <- experiment_site_keys(sites)
    ranges
}

#' Identify reference sites including their measurement geometry
#'
#' @param sites Data frame of genomic anchors and measurement annotations.
#' @return Character site IDs.
#' @noRd
experiment_site_keys <- function(sites) {
    paste(
        sites$seqname, 
        sites$position, 
        sites$strand, 
        sites$context,
        sites$modification, 
        sep = ":", 
        recycle0 = TRUE
    )
}

#' Populate the matrices for the ranged experiment
#'
#' @param measurements Combined reference calls.
#' @param sites Named `GRanges` for the experiment rows.
#' @param samples Experiment column names.
#' @return A named list of assay matrices. Duplicate cells signal an error.
#' @noRd
experiment_assays <- function(measurements, sites, samples) {
    rows <- match(experiment_site_keys(measurements), names(sites))
    columns <- match(measurements$sample, samples)
    duplicate <- anyDuplicated(data.frame(rows, columns))
    if (duplicate) {
        site <- names(sites)[rows[[duplicate]]]
        sample <- samples[columns[[duplicate]]]
        cli::cli_abort(
            paste(
                "Multiple reference calls for sample {.val {sample}}",
                "at {.val {site}}."
            ),
            class = "methylVCF_duplicate_site",
            sample = sample, site = site, call = NULL
        )
    }
    cells <- cbind(rows, columns)
    columns_to_keep <- c("fraction", "informative_depth", "modified_depth", "coverage_status")
    values <- measurements[columns_to_keep]
    names(values)[[1L]] <- "reported_fraction"
    observed <- which(values$coverage_status == "observed")
    fraction <- rep(NA_real_, nrow(measurements))
    fraction[observed] <- values$reported_fraction[observed]
    lapply(c(list(fraction = fraction), values), function(value) {
        assay <- matrix(
            value[NA_integer_], 
            nrow = length(sites), 
            ncol = length(samples)
        )
        assay[cells] <- value
        assay
    })
}

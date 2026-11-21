#' Query reported 5mC in a genomic region
#'
#' Interpret 5mC fields from DRAGEN and Rastair and retain evidence supplied
#' by the file, including explicit missing slots and zero fractions. The
#' query filters observations by their reported anchors after decoding.
#'
#' @param file Path to a local text VCF or gVCF, optionally compressed with
#'   gzip or BGZF.
#' @param which A `GRanges` with 1-based inclusive regions to query. The
#'   reader combines overlapping regions and ignores strand in the query.
#' @param format The producer's format, `"dragen"` (default) or `"rastair"`.
#'   Use the `format` reported by [inspectMethylVcf()] when identified.
#' @param samples Character vector of VCF sample names to return, or `NULL`
#'   to retain all samples.
#' @param on_problem `"error"` stops at the first unsupported or invalid
#'   record, sample, or measurement. `"report"` retains readable observations
#'   and records omissions in the result metadata. Header, file, and index
#'   errors always stop the import.
#'
#' @return A `GRanges` with one row per reported measurement slot and a
#'   one-base range at its reported anchor. Metadata columns include the
#'   sample, context, scope, fraction, informative and modified read depths,
#'   genotype, optional haplotype and phase-set identifiers, allele index
#'   (0 for REF and 1 for the first ALT), and the source record coordinates.
#'   `allele_offset` is the zero-based offset within the allele sequence or
#'   reference block. In a gVCF block, `allele` gives the sequence at the
#'   record's anchor only. `methylation_mask` retains the original DRAGEN
#'   `INFO/M5mC` string in REF/ALT order when present.
#'   The `Seqinfo` retains all contigs, lengths and assembly tags declared in
#'   the VCF header. Queried contigs absent from the header have unknown
#'   lengths and assemblies.
#'
#' @details DRAGEN support includes reference blocks with `INFO/END`, diploid
#'   SNVs including multiallelic sites, and measurements on mappable alleles
#'   of complex records. When REF and ALT differ in length, projection is
#'   limited to their shared leading anchor. Rastair reference and biallelic
#'   single-base contexts use their producer-specific order. When both
#'   `CPG` and `CPGnovo` are set, values must be separable by context.
#'   Rastair 2.2.0 can emit a single slot for such records. The reader
#'   reports that representation as unsupported.
#'   A Rastair slot with more modified reads (`ADM5mC`) than informative
#'   reads (`DPM5mC`) is invalid when both depths are available.
#'
#'   `scope` is `"pooled"`, `"allele"`, `"haplotype"`, or `"context"`. A DRAGEN
#'   phased homozygote with one measurement per site has pooled scope and a
#'   missing haplotype. A Rastair allele assignment identifies a local
#'   sequence allele to which the CpG context can be attributed. Its haplotype
#'   is missing. If multiple local alleles could anchor the context, the
#'   reader retains the measurement with a missing allele index and
#'   unspecified strand. Unreported measurements remain missing. The file-level
#'   metadata include the chosen format, selected sample IDs (even for
#'   samples without calls), and available header provenance.
#'
#'   `coverage_status` is `"observed"` for a reported fraction with positive
#'   informative depth. `"no_informative_depth"` takes precedence when the
#'   depth is zero. `"missing_fraction"` indicates a missing fraction and
#'   `"unknown_depth"` indicates a reported fraction with missing depth.
#'   `fraction` retains the producer's value.
#'
#'   The reader uses Rsamtools for BGZF files with a companion TBI or CSI
#'   index, preferring CSI when both are present. Otherwise it scans the
#'   input sequentially in chunks. Indexed queries retrieve the requested
#'   intervals and decode a source record only in its first overlapping
#'   interval, keeping all measurements within the query.
#'   Anchors include the start of a collapsed CpG dyad. Retained rows preserve
#'   source order within each record.
#'
#' @seealso [inspectMethylVcf()] for checking the header and index,
#'   [asMethylExperiment()] for comparing reference sites across samples,
#'   and `vignette("methylVCF")` for the complete workflow.
#'
#' @section Import diagnostics:
#'   `S4Vectors::metadata(result)$diagnostics` contains one row per import
#'   issue, with `chrom`, `pos`, `sample` when known, `type` (`"invalid"` or
#'   `"unsupported"`), `allele_index`, `allele_offset`, and `value_index`
#'   when available, and `reason`. In `"error"` mode, the first such issue
#'   signals a `methylVCF_import_error` condition with those fields. In
#'   `"report"` mode, a single `methylVCF_import_warning` condition carries the
#'   total issue count in `n_issues`.
#'
#' @examples
#' file <- system.file("extdata", "dragen-example.g.vcf", package = "methylVCF")
#' region <- GenomicRanges::GRanges("chr21", IRanges::IRanges(13204474L))
#' result <- readMethylVcf(file, region, format = "dragen")
#' result
#'
#' @testexamples
#' testthat::expect_s4_class(result, "GRanges")
#' measurements <- as.data.frame(result)
#' testthat::expect_identical(as.character(measurements$seqnames), "chr21")
#' testthat::expect_identical(measurements$start, 13204474L)
#' testthat::expect_identical(measurements$record_start, 13204471L)
#' testthat::expect_equal(measurements$fraction, 0.2381)
#' testthat::expect_identical(measurements$informative_depth, 21L)
#'
#' @export
#' @md
readMethylVcf <- function(
    file,
    which,
    format = c("dragen", "rastair"),
    samples = NULL,
    on_problem = c("error", "report")
) {
    format <- match.arg(format)
    on_problem <- match.arg(on_problem)
    validate_query_arguments(file, which, samples)

    regions <- GenomicRanges::reduce(which, ignore.strand = TRUE)
    coordinates <- as.data.frame(regions)[c("seqnames", "start", "end")]
    by_contig <- split(coordinates, coordinates$seqnames, drop = TRUE)
    file_name <- basename(file)
    index <- if (length(regions)) find_vcf_index(file)
    access <- if (is.null(index)) "sequential" else "tabix"

    con <- open_vcf_text(file)

    on.exit(close(con), add = TRUE)
    header <- read_vcf_header(con)
    selected <- samples %||% header$samples
    if (!all(selected %in% header$samples)) {
        absent <- paste(setdiff(selected, header$samples), collapse = ", ")
        cli::cli_abort("Samples absent from VCF: {absent}", call = NULL)
    }

    decode <- switch(format, dragen = decode_dragen, rastair = decode_rastair)
    reporter <- new_problem_reporter(on_problem)
    calls <- read_call_batches(
        con, by_contig, header, selected, decode, reporter$add,
        file = file, index = index
    )
    assemble_methyl_ranges(
        calls,
        coordinates$seqnames,
        header,
        reporter$diagnostics(),
        format,
        file_name,
        access,
        selected
    )
}

#' Validate the local input and select query arguments
#'
#' @param file Path to the local input file.
#' @param which Requested genomic ranges.
#' @param samples Requested sample names or `NULL`.
#' @return Invisibly `NULL` after validating the supplied arguments.
#' @noRd
validate_query_arguments <- function(file, which, samples) {
    validate_vcf_file(file)
    valid_ranges <- inherits(which, "GRanges") && all(IRanges::start(which) >= 1L & IRanges::width(which) >= 1L)
    if (!valid_ranges) {
        cli::cli_abort(
            "'which' must be a GRanges with positive, nonempty coordinates.",
            call = NULL
        )
    }
    if (!is.null(samples) && !is_valid_sample_ids(samples)) {
        cli::cli_abort(
            "'samples' must contain unique, nonempty sample names.",
            call = NULL
        )
    }
    invisible(NULL)
}

#' Validate a local VCF path
#' @param file Path to a local VCF or gVCF.
#' @return Invisibly `NULL` for a local file path.
#' @noRd
validate_vcf_file <- function(file) {
    valid_file <- is.character(file) && length(file) == 1L && !is.na(file) && file.exists(file) && !dir.exists(file)
    if (!valid_file) {
        cli::cli_abort("'file' must be a path to a local VCF file.", call = NULL)
    }
    invisible(NULL)
}

#' Open a text connection to a VCF
#' @param file Path to a plain, gzip or BGZF VCF.
#' @return An open connection positioned at the start of the file.
#' @noRd
open_vcf_text <- function(file) {
    if (grepl("[.](gz|bgz)$", file, ignore.case = TRUE)) { return(gzfile(file, "rt")) }
    base::file(file, "rt")
}

#' Decode an open VCF connection in batches
#'
#' @param con Text connection positioned after the header.
#' @param by_contig Reduced query intervals keyed by contig.
#' @param header Parsed VCF header.
#' @param selected Sample names requested by the caller.
#' @param decode Producer-specific decoder for one record and sample.
#' @param problem Callback for invalid or unsupported measurements.
#' @param file Source VCF path.
#' @param index TBI or CSI path for indexed access, or `NULL` for sequential.
#' @return A combined table with the schema from `empty_result()`.
#' @noRd
read_call_batches <- function(
    con,
    by_contig,
    header,
    selected,
    decode,
    problem,
    file,
    index
) {
    batches <- list(empty_result())
    if (!length(by_contig) || !length(selected)) { return(batches[[1L]]) }
    sample_columns <- 9L + match(selected, header$samples)
    names(sample_columns) <- names(selected)
    column_count <- 9L + length(header$samples)
    decode_batch <- function(lines, previous_end = 0L) {
        rows <- lapply(
            lines[nzchar(lines)],
            decode_line,
            by_contig = by_contig,
            samples = selected,
            sample_columns = sample_columns,
            column_count = column_count,
            decode = decode,
            problem = problem,
            previous_end = previous_end
        )
        do.call(rbind, rows)
    }
    if (!is.null(index)) {
        indexed <- read_indexed_batches(file, index, by_contig, decode_batch)
        return(do.call(rbind, c(batches, indexed)))
    }
    repeat {
        lines <- readLines(con, n = 4096L, warn = FALSE)
        if (!length(lines)) { break }
        batches[[length(batches) + 1L]] <- decode_batch(lines)
    }
    do.call(rbind, batches)
}

#' Assemble the public ranges and their provenance
#'
#' @param calls Combined methylation call table.
#' @param requested_contigs Contigs from the requested query.
#' @param header Parsed VCF header.
#' @param diagnostics Diagnostic rows for this import.
#' @param format Producer name used for decoding.
#' @param file_name Basename of the requested input file.
#' @param access Indexed or sequential access method.
#' @param samples VCF sample IDs selected for this query.
#' @return A `GRanges` with provenance and diagnostics. If any issues were
#'   recorded, it emits one structured warning.
#' @noRd
assemble_methyl_ranges <- function(
    calls,
    requested_contigs,
    header,
    diagnostics,
    format,
    file_name,
    access,
    samples
) {
    contigs <- unique(c(as.character(requested_contigs), calls$chrom))
    sequence_info <- header$seqinfo
    Seqinfo::seqlevels(sequence_info) <- union(Seqinfo::seqlevels(sequence_info), contigs)
    out <- GenomicRanges::makeGRangesFromDataFrame(
        calls,
        keep.extra.columns = TRUE,
        seqnames.field = "chrom",
        start.field = "position",
        end.field = "position",
        strand.field = "strand",
        seqinfo = sequence_info
    )

    S4Vectors::metadata(out) <- list(
        format = format,
        source = header$source,
        caller_version = header$caller_version,
        vcf_version = header$version,
        file = file_name,
        reference = header$reference,
        access = access,
        samples = samples,
        coordinate = "reported base anchor",
        diagnostics = diagnostics
    )

    n_issues <- nrow(diagnostics)
    if (n_issues) {
        cli::cli_warn(
            paste(
                "The import reported {n_issues} issue{?s}.",
                "See metadata(result)$diagnostics."
            ),
            class = "methylVCF_import_warning", n_issues = n_issues, call = NULL
        )
    }
    out
}

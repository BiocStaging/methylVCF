#' VCF FORMAT fields carrying methylation fractions and read depths
#' @noRd
methyl_format_fields <- c("M5mC", "DPM5mC", "ADM5mC")

#' Read one VCF field with the missing-value default
#'
#' @param fields Named list of INFO or sample FORMAT values.
#' @param name Field name to look up.
#' @return The stored value, including `TRUE` for an INFO flag, or `"."`
#'   for an absent or empty field.
#' @noRd
field <- function(fields, name) {
    value <- fields[[name]]
    if (is.null(value) || !nzchar(value)) "." else value
}

#' Check whether sample identifiers can be used as unique keys
#'
#' @param ids Candidate sample identifiers.
#' @return A single logical value. An empty character vector is valid.
#' @noRd
is_valid_sample_ids <- function(ids) {
    is.character(ids) && !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids)
}

#' Validate an exact VCF FILTER selection
#'
#' @param filters Character vector of complete FILTER values, or `NULL` to
#'   retain all values.
#' @return Invisibly `NULL` for a valid selection.
#' @noRd
validate_methyl_filters <- function(filters) {
    filters <- filters %||% character()
    valid <- is.character(filters) && !anyNA(filters) && all(nzchar(filters))
    if (!valid) {
        cli::cli_abort(
            "{.arg filters} must contain VCF FILTER values or be NULL.",
            call = NULL
        )
    }
    invisible(NULL)
}

#' Split comma-separated VCF values
#'
#' @param value A single character string containing the field value.
#' @return A character vector in source order. A trailing empty element is
#'   retained so later validation can distinguish it from a missing `"."`.
#' @noRd
split_values <- function(value) {
    # Append a sentinel to preserve a malformed trailing empty element.
    parts <- strsplit(paste0(value, ",."), ",", fixed = TRUE)[[1L]]
    parts[-length(parts)]
}

#' Parse methylation fractions and depths from a sample
#'
#' @param fields Named list of sample FORMAT values.
#' @return A list of numeric `fraction` and integer `informative_depth` and
#'   `modified_depth` vectors in source order. Dots become typed `NA`s.
#'   Fractions must lie in [0, 1], and depths must be nonnegative integers.
#' @noRd
read_methyl_values <- function(fields) {
    parse_value <- function(name, is_depth = FALSE) {
        value <- split_values(field(fields, name))
        missing <- value == "."
        parsed <- utils::type.convert(value, na.strings = ".", as.is = TRUE)
        if (!is.numeric(parsed) && !all(missing)) {
            abort_methyl_record(paste0("Invalid ", name, " value."))
        }
        parsed <- as.double(parsed)
        invalid_number <- any(!missing & (!is.finite(parsed) | parsed < 0))
        invalid_fraction <- !is_depth && any(parsed > 1, na.rm = TRUE)
        if (invalid_number || invalid_fraction) {
            abort_methyl_record(paste0("Invalid ", name, " value."))
        }
        if (!is_depth) { return(parsed) }

        invalid_depth <- any(!missing & !grepl("^[0-9]+$", value)) || any(parsed > .Machine$integer.max, na.rm = TRUE)
        if (invalid_depth) {
            abort_methyl_record(paste0("Invalid ", name, " read depth."))
        }
        as.integer(parsed)
    }

    list(
        fraction = parse_value("M5mC"),
        informative_depth = parse_value("DPM5mC", is_depth = TRUE),
        modified_depth = parse_value("ADM5mC", is_depth = TRUE)
    )
}

#' Align a parsed FORMAT vector to the caller's measurement slots
#'
#' @param values Fraction or depth values for one sample field.
#' @param n Number of slots reported by the caller.
#' @param reason Caller-specific diagnostic for an incompatible count.
#' @return Values in slot order. A wholly unreported field expands to typed
#'   missing values for each slot.
#' @noRd
align_methyl_values <- function(values, n, reason) {
    if (!length(values) || (length(values) == 1L && is.na(values))) {
        return(rep(values[NA_integer_], n))
    }
    if (length(values) != n) { abort_methyl_record(reason) }
    values
}

#' Construct the shared methylation call table
#'
#' @param position Integer vector of 1-based genomic anchors. Its length
#'   determines the number of rows and the length of the other columns.
#' @param strand Strand symbol, with `"*"` for a combined or unknown strand.
#' @param context Sequence context, such as `"CG"` or `"CHH"`.
#' @param scope Evidence level for the measurement.
#' @param fraction Reported 5mC fraction on the [0, 1] scale.
#' @param informative_depth,modified_depth Reported read counts, with `NA`
#'   for unreported counts.
#' @param allele_index Zero-based index into REF followed by ALT alleles.
#' @param haplotype One-based genotype slot for a phased measurement.
#' @param motif_origin Producer-specific origin of the sequence context.
#' @param value_index One-based index of the measurement in the FORMAT field.
#' @param allele_offset Zero-based offset within the allele or reference block.
#' @return A typed data frame with one row per anchor. Columns are recycled
#'   to the number of anchors. `coverage_status` distinguishes a reported
#'   fraction with positive depth, depth zero, and missing values.
#' @noRd
new_calls <- function(
    position,
    strand,
    context,
    scope,
    fraction,
    informative_depth,
    modified_depth = NA_integer_,
    allele_index = NA_integer_,
    haplotype = NA_integer_,
    motif_origin = NA_character_,
    value_index = seq_along(position),
    allele_offset = 0L
) {
    n <- length(position)
    calls <- data.frame(
        position = as.integer(position),
        strand = rep_len(as.character(strand), n),
        context = rep_len(as.character(context), n),
        modification = rep("5mC", n),
        scope = rep_len(as.character(scope), n),
        fraction = rep_len(as.numeric(fraction), n),
        informative_depth = rep_len(as.integer(informative_depth), n),
        modified_depth = rep_len(as.integer(modified_depth), n),
        allele_index = rep_len(as.integer(allele_index), n),
        haplotype = rep_len(as.integer(haplotype), n),
        motif_origin = rep_len(as.character(motif_origin), n),
        value_index = rep_len(as.integer(value_index), n),
        allele_offset = rep_len(as.integer(allele_offset), n),
        stringsAsFactors = FALSE
    )

    calls$coverage_status <- rep("observed", n)
    calls$coverage_status[is.na(calls$informative_depth)] <- "unknown_depth"
    calls$coverage_status[is.na(calls$fraction)] <- "missing_fraction"
    no_depth <- which(calls$informative_depth == 0L)
    calls$coverage_status[no_depth] <- "no_informative_depth"
    calls
}

#' Construct an empty decoder result
#'
#' @return A zero-row data frame with the column types from `new_calls()`.
#' @noRd
empty_calls <- function() {
    new_calls(
        integer(), character(), character(), character(), numeric(), integer()
    )
}

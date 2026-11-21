#' Decode a VCF record and retain its queried measurement anchors
#'
#' @param line A single tab-delimited data record.
#' @param by_contig Reduced query intervals keyed by contig, with sorted
#'   1-based inclusive starts and ends.
#' @param samples Selected sample IDs, in their requested order.
#' @param sample_columns VCF column numbers in the same order as `samples`.
#' @param column_count Number of columns expected from the original header.
#' @param decode Decoder called with the parsed record and one sample's fields.
#' @param problem Callback that aborts in strict mode or records an issue
#'   and returns `NULL` in report mode.
#' @param previous_end End of the previous indexed query window, or zero.
#' @return A call table with source and sample metadata, or `NULL` if there
#'   is no overlapping record or no selected sample supplies a call table.
#' @noRd
decode_line <- function(
    line,
    by_contig,
    samples,
    sample_columns,
    column_count,
    decode,
    problem,
    previous_end = 0L
) {
    cols <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    chrom <- cols[[1L]]
    ranges <- by_contig[[chrom]]
    if (is.null(ranges)) { return(NULL) }

    pos <- NA_integer_
    record <- tryCatch(
        {
            if (length(cols) >= 2L) { pos <- as_position(cols[[2L]]) }
            # An indexed hit reaching this window also spans the previous
            # window if it starts at or before the previous window's end.
            if (!is.na(pos) && pos <= previous_end) { return(NULL) }
            parse_vcf_record(cols, chrom, pos, ranges, column_count)
        },
        methylVCF_record_error = function(cnd) {
            problem(chrom, pos, NA_character_, "invalid", cnd$reason)
        }
    )
    if (is.null(record)) { return(NULL) }

    rows <- Map(
        function(column, sample) {
            decode_sample_calls(
                cols[[column]], sample, record, ranges, decode, problem
            )
        },
        sample_columns,
        samples
    )
    do.call(rbind, rows)
}

#' Parse one source record within its query intervals
#'
#' @param cols Tab-delimited fields from one VCF line.
#' @param chrom Contig name.
#' @param pos Validated POS or `NA_integer_` if the column is missing.
#' @param ranges Query intervals for the record's contig.
#' @param column_count Number of columns expected from the header.
#' @return A parsed record or `NULL` for a nonoverlapping record. Malformed
#'   records signal an error for the caller to report.
#' @noRd
parse_vcf_record <- function(cols, chrom, pos, ranges, column_count) {
    if (!is.na(pos) && pos > ranges$end[[nrow(ranges)]]) { return(NULL) }
    if (length(cols) != column_count) {
        abort_methyl_record("Record column count differs from VCF header.")
    }

    info <- parse_info(cols[[8L]])
    ref <- cols[[4L]]
    end <- vcf_record_end(pos, ref, info$END)
    if (!any(end >= ranges$start & pos <= ranges$end)) { return(NULL) }

    alt <- if (cols[[5L]] == ".") character() else split_values(cols[[5L]])
    keys <- strsplit(cols[[9L]], ":", fixed = TRUE)[[1L]]
    if (anyDuplicated(keys)) { abort_methyl_record("Duplicate FORMAT keys.") }

    list(
        chrom = chrom,
        pos = pos,
        end = end,
        ref = ref,
        alt = alt,
        info = info,
        filter = cols[[7L]],
        keys = keys
    )
}

#' Determine a VCF reference span without overflowing an integer
#'
#' @param pos Validated 1-based VCF POS.
#' @param ref REF sequence supplied by the record.
#' @param end Explicit INFO/END string or `NULL` for an implicit span.
#' @return Validated 1-based end as an integer.
#' @noRd
vcf_record_end <- function(pos, ref, end) {
    if (!nzchar(ref)) { abort_methyl_record("Empty VCF REF.") }
    ref_end <- as.double(pos) + nchar(ref) - 1
    if (ref_end > .Machine$integer.max) {
        abort_methyl_record(
            "Record END exceeds the supported coordinate range."
        )
    }
    if (is.null(end)) { return(as.integer(ref_end)) }
    end <- as_position(end)
    if (end < pos) { abort_methyl_record("Record END precedes POS.") }
    if (end < ref_end) { abort_methyl_record("Record END precedes the last REF base.") }
    end
}

#' Decode one sample in an overlapping source record
#'
#' @param column Colon-delimited sample column from the VCF record.
#' @param sample Header name of that sample.
#' @param record Parsed record, including the ordered FORMAT keys.
#' @param ranges Reduced query intervals for the record's contig.
#' @param decode Producer-specific decoder.
#' @param problem Callback to report a sample or measurement issue.
#' @return Calls anchored within the requested intervals, with source and
#'   sample metadata, or `NULL` if the decoder returns no calls.
#' @noRd
decode_sample_calls <- function(column, sample, record, ranges, decode, problem) {
    decoded <- tryCatch(
        {
            values <- strsplit(column, ":", fixed = TRUE)[[1L]]
            keys <- record$keys
            if (length(values) > length(keys)) {
                abort_methyl_record("More sample values than FORMAT keys.")
            }
            fields <- stats::setNames(as.list(values), keys[seq_along(values)])
            list(calls = decode(record, fields), fields = fields)
        },
        methylVCF_record_error = function(cnd) {
            unsupported <- inherits(cnd, "methylVCF_unsupported")
            type <- if (unsupported) "unsupported" else "invalid"
            problem(record$chrom, record$pos, sample, type, cnd$reason)
        }
    )

    calls <- decoded$calls
    skipped <- attr(calls, "issues")
    if (!is.null(skipped)) {
        for (index in seq_len(nrow(skipped))) {
            issue <- skipped[index, ]
            problem(
                record$chrom,
                record$pos,
                sample,
                issue$type,
                issue$reason,
                issue$allele_index,
                issue$allele_offset,
                issue$value_index
            )
        }
    }
    if (is.null(calls) || !nrow(calls)) { return(NULL) }

    hits <- findInterval(calls$position, ranges$start)
    keep <- calls$position <= c(0L, ranges$end)[hits + 1L]
    attach_record(calls[keep, , drop = FALSE], record, decoded$fields, sample)
}

#' Parse a 1-based VCF coordinate
#'
#' @param value Coordinate as a string of decimal digits.
#' @return A positive integer. Unrepresentable or malformed values signal
#'   a record error.
#' @noRd
as_position <- function(value) {
    scalar <- is.character(value) && length(value) == 1L && !is.na(value)
    if (!scalar || !grepl("^[0-9]+$", value)) {
        abort_methyl_record("Invalid VCF coordinate.")
    }
    parsed <- as.double(value)
    if (parsed < 1 || parsed > .Machine$integer.max) {
        abort_methyl_record("Invalid VCF coordinate.")
    }
    as.integer(parsed)
}

#' Parse VCF INFO values and flags
#'
#' @param info A semicolon-delimited INFO field, or `"."` if it is missing.
#' @return A list keyed by INFO name. Assignments keep their raw strings,
#'   bare flags become `TRUE`, and missing INFO gives an empty list.
#' @noRd
parse_info <- function(info) {
    if (info == ".") { return(list()) }
    parts <- strsplit(info, ";", fixed = TRUE)[[1L]]
    keys <- sub("=.*$", "", parts)
    if (anyDuplicated(keys)) { abort_methyl_record("Duplicate INFO keys.") }
    values <- lapply(parts, function(part) {
        if (grepl("=", part, fixed = TRUE)) sub("^[^=]*=", "", part) else TRUE
    })
    stats::setNames(values, keys)
}

#' Attach sample and source record metadata to decoded calls
#'
#' A gVCF block's REF describes its anchor base. Calls inside the block keep
#' their allele index and offset and have a missing allele sequence.
#'
#' @param calls Decoder call table with genomic anchors and allele indices.
#' @param record Parsed source record, including boundaries, alleles and filter.
#' @param fields Named list of the selected sample's FORMAT fields.
#' @param sample Name of the sample in the VCF header.
#' @return The call table with sample, genotype, phase set and source columns.
#' @noRd
attach_record <- function(calls, record, fields, sample) {
    n <- nrow(calls)
    calls$chrom <- rep(record$chrom, n)
    calls$sample <- rep(sample, n)
    calls$genotype <- rep(field(fields, "GT"), n)
    calls$genotype[calls$genotype == "."] <- NA_character_
    calls$phase_set <- rep(field(fields, "PS"), n)
    calls$phase_set[calls$phase_set == "."] <- NA_character_

    calls$allele <- c(record$ref, record$alt)[calls$allele_index + 1L]
    if (!is.null(record$info$END) && nchar(record$ref) == 1L) {
        calls$allele[calls$position != record$pos] <- NA_character_
    }

    calls$record_start <- rep(record$pos, n)
    calls$record_end <- rep(record$end, n)
    calls$ref <- rep(record$ref, n)
    calls$alt <- rep(paste(record$alt, collapse = ","), n)
    calls$filter <- rep(record$filter, n)
    mask <- record$info[["M5mC"]]
    mask <- if (is.character(mask)) mask else NA_character_
    calls$methylation_mask <- rep(mask, n)
    calls
}

#' Construct an empty, fully annotated reader result
#'
#' @return A zero-row data frame with the call and source metadata columns
#'   used to combine decoded batches and construct the final `GRanges`.
#' @noRd
empty_result <- function() {
    record <- list(
        chrom = "",
        pos = 1L,
        end = 1L,
        ref = "N",
        alt = character(),
        filter = "."
    )
    attach_record(empty_calls(), record, list(), "")
}

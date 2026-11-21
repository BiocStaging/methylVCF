#' Audit REF and ALT methylation at the same DRAGEN anchor
#'
#' @param calls Regional `GRanges` from `readMethylVcf()` in DRAGEN format,
#'   including the original INFO/M5mC masks. Use one import at a time.
#' @param filters Exact VCF FILTER values to retain, or `NULL` for all records.
#' @return A `GRanges` with one row per observed anchor, sample and source
#'   record for heterozygous biallelic substitutions. The range is unstranded
#'   and the two allele-specific columns retain the original strands, depths,
#'   fractions, coverage statuses and local haplotype slots. `pair_status`
#'   records why a comparison can or cannot be made. `delta` is ALT minus REF
#'   for comparable pairs and `NA_real_` otherwise. Stored masks distinguish
#'   a caller-unreported counterpart from a missing decoded row.
#'
#' @details Real SNVs and equal-length substitutions are eligible for the
#'   audit. A delta requires the same explicitly preserved CpG sequence and
#'   measurement geometry on both alleles. The two bases of the CpG must be
#'   present in both source sequences. Phase slots are local to their record
#'   and phase set. A dot mask can mean unreported or not applicable. It
#'   cannot establish the biological loss of a CpG or zero methylation.
#'
#'   Status values distinguish a comparable pair, insufficient coverage at a
#'   shared CpG, different measurement geometry, an unverified CpG, one allele
#'   reported (`ref_only` or `alt_only`), an absent expected call, and
#'   ambiguous duplicates.
#' @noRd
pairAllelicMethylation <- function(calls, filters = c("PASS", ".")) {
    validate_methyl_filters(filters)
    dragen <- inherits(calls, "GRanges") && identical(S4Vectors::metadata(calls)$format, "dragen")
    if (!dragen) {
        cli::cli_abort(
            "{.arg calls} must be GRanges from a DRAGEN import.", call = NULL
        )
    }
    fields <- S4Vectors::mcols(calls)
    record_fields <- c(
        "sample", "record_start", "record_end", "ref", "alt", "genotype",
        "phase_set", "filter", "methylation_mask", "modification",
        "allele_offset"
    )
    slot_fields <- c(
        "context", "fraction", "informative_depth", "modified_depth",
        "coverage_status", "haplotype"
    )
    required <- c(record_fields, slot_fields, "allele_index", "scope")
    if (!all(required %in% names(fields))) {
        cli::cli_abort(
            paste(
                "Read {.arg calls} again with readMethylVcf()",
                "to retain the source masks."
            ),
            call = NULL
        )
    }
    eligible <- fields$scope %in% c("allele", "haplotype")
    if (!is.null(filters)) { eligible <- eligible & fields$filter %in% filters }
    calls <- calls[eligible]
    fields <- fields[eligible, , drop = FALSE]
    variant <- allelic_substitutions(fields)
    calls <- calls[variant$rows]
    fields <- fields[variant$rows, , drop = FALSE]
    expected_anchor <- as.double(fields$record_start) + fields$allele_offset
    original_anchor <- GenomicRanges::start(calls) == expected_anchor
    if (!isTRUE(all(GenomicRanges::width(calls) == 1L & original_anchor))) {
        cli::cli_abort(
            "Candidate ranges must retain their original one-base anchors.",
            call = NULL
        )
    }

    record <- fields[, record_fields, drop = FALSE]
    site <- S4Vectors::DataFrame(
        seqname = as.character(GenomicRanges::seqnames(calls)),
        anchor = GenomicRanges::start(calls), record
    )
    groups <- S4Vectors::selfmatch(site)
    first <- which(!duplicated(groups))
    ref_rows <- which(fields$allele_index == 0L)
    alt_rows <- which(fields$allele_index == variant$alt_index)
    ref <- allelic_slots(
        calls, fields, slot_fields, ref_rows, first, groups, prefix = "ref"
    )
    alt <- allelic_slots(
        calls, fields, slot_fields, alt_rows, first, groups, prefix = "alt"
    )

    info <- record[first, , drop = FALSE]
    info$record_alt <- info$alt
    info$alt <- variant$alt[first]
    info$alt_index <- variant$alt_index[first]
    masks <- allelic_mask_at_offset(info)
    paired <- cbind(info, masks, ref$slots, alt$slots)
    paired$pair_status <- allelic_pair_status(paired, ref$count, alt$count)
    comparable <- which(paired$pair_status == "comparable")
    paired$delta <- rep(NA_real_, nrow(paired))
    paired$delta[comparable] <- paired$alt_fraction[comparable] - paired$ref_fraction[comparable]

    result <- calls[first]
    GenomicRanges::strand(result) <- rep("*", length(result))
    S4Vectors::mcols(result) <- paired
    S4Vectors::metadata(result)["pair_filters"] <- list(filters)
    result
}

#' Identify called biallelic, equal-length substitutions
#'
#' @param fields Metadata columns for candidate DRAGEN calls.
#' @return Row indices and the original ALT indices and sequences. Only
#'   gVCF placeholders are ignored when testing for biallelic records.
#' @noRd
allelic_substitutions <- function(fields) {
    alternates <- unique(fields$alt)
    alleles <- strsplit(alternates, ",", fixed = TRUE)
    alt_index <- vapply(alleles, function(alleles) {
        literal <- which(!alleles %in% c("<NON_REF>", "<*>"))
        if (length(literal) == 1L) literal else NA_integer_
    }, integer(1L))
    sequence <- vapply(seq_along(alleles), function(index) {
        alleles[[index]][alt_index[index]]
    }, character(1L))

    alt <- toupper(sequence[match(fields$alt, alternates)])
    ref <- toupper(fields$ref)
    index <- alt_index[match(fields$alt, alternates)]
    genotypes <- cbind(
        paste0("0/", index, recycle0 = TRUE),
        paste0(index, "/0", recycle0 = TRUE),
        paste0("0|", index, recycle0 = TRUE),
        paste0(index, "|0", recycle0 = TRUE)
    )
    heterozygous <- rowSums(genotypes == fields$genotype, na.rm = TRUE) > 0L
    literal <- grepl("^[ACGT]+$", ref) & grepl("^[ACGT]+$", alt)
    substitution <- ref != alt & nchar(ref) == nchar(alt)
    ref_span <- fields$record_end - fields$record_start == nchar(ref) - 1L
    called <- fields$allele_index == 0L | fields$allele_index == index
    eligible <- !is.na(index) & heterozygous & literal & substitution
    rows <- which(eligible & ref_span & called)
    list(rows = rows, alt = alt[rows], alt_index = index[rows])
}

#' Align one allele's slots to the first row of each source group
#'
#' @param calls Candidate ranges.
#' @param fields Corresponding observation metadata.
#' @param slot_fields Columns to retain for each allele.
#' @param rows Indices for the requested allele.
#' @param first First row for each source group.
#' @param groups First-row index for each candidate range.
#' @param prefix Allele prefix for the returned columns.
#' @return Slot columns and counts per group. Absent and ambiguous slots
#'   get typed missing values.
#' @noRd
allelic_slots <- function(
    calls,
    fields,
    slot_fields,
    rows,
    first,
    groups,
    prefix
) {
    count <- tabulate(match(groups[rows], first), nbins = length(first))
    slot_rows <- rows[match(first, groups[rows])]
    slot_rows[count != 1L] <- NA_integer_
    strand <- as.character(GenomicRanges::strand(calls))[slot_rows]
    slots <- S4Vectors::DataFrame(
        strand = strand, 
        fields[slot_rows, slot_fields, drop = FALSE]
    )
    names(slots) <- paste0(prefix, "_", names(slots))
    list(slots = slots, count = count)
}

#' Recover the two DRAGEN mask characters at a paired anchor
#'
#' @param info Source records including the full mask, offset and ALT index.
#' @return A `DataFrame` with one mask character per allele, or `NA` when
#'   the character cannot be recovered.
#' @noRd
allelic_mask_at_offset <- function(info) {
    masks <- strsplit(info$methylation_mask, ",", fixed = TRUE)
    code <- function(indices) {
        vapply(seq_along(masks), function(row) {
            allele <- masks[[row]][indices[row]]
            position <- info$allele_offset[row] + 1L
            if (identical(allele, ".")) { return(".") }
            value <- substring(allele, position, position)
            if (value %in% c(".", "Z", "z", "x", "h")) value else NA_character_
        }, character(1L))
    }
    S4Vectors::DataFrame(
        ref_mask = code(rep(1L, nrow(info))),
        alt_mask = code(info$alt_index + 1L)
    )
}

#' Classify candidate comparisons without imputing missing evidence
#'
#' @param pairs Paired source and allele columns.
#' @param ref_count,alt_count Number of source observations for each allele.
#' @return A status string per candidate. Geometry and coverage are evaluated
#'   for unique, present allele slots.
#' @noRd
allelic_pair_status <- function(pairs, ref_count, alt_count) {
    status <- rep("missing_call", nrow(pairs))
    ref_only <- ref_count == 1L & alt_count == 0L & pairs$alt_mask == "."
    alt_only <- alt_count == 1L & ref_count == 0L & pairs$ref_mask == "."
    status[which(ref_only)] <- "ref_only"
    status[which(alt_only)] <- "alt_only"
    status[ref_count > 1L | alt_count > 1L] <- "ambiguous"
    paired <- which(ref_count == 1L & alt_count == 1L)
    both <- pairs[paired, , drop = FALSE]

    same_geometry <- both$ref_context == both$alt_context & both$ref_strand == both$alt_strand & both$ref_mask == both$alt_mask
    status[paired] <- "different_geometry"
    aligned <- which(same_geometry)
    status[paired[aligned]] <- "unverified_cpg"
    verified <- intersect(aligned, which(has_shared_cpg(both)))
    status[paired[verified]] <- "insufficient_coverage"
    observed <- both$ref_coverage_status == "observed" & both$alt_coverage_status == "observed"
    finite <- is.finite(both$ref_fraction) & is.finite(both$alt_fraction)
    covered <- both$ref_informative_depth > 0L & both$alt_informative_depth > 0L
    comparable <- intersect(verified, which(observed & finite & covered))
    status[paired[comparable]] <- "comparable"
    status
}

#' Verify the CpG sequence and reported geometry of two aligned allele slots
#'
#' @param pairs Source records with a unique REF slot and a unique ALT slot.
#' @return Logical vector identifying the same CpG in both literal sequences.
#'   The caller checks agreement between the two strands, masks and contexts.
#' @noRd
has_shared_cpg <- function(pairs) {
    dyad <- pairs$ref_mask == "Z" & pairs$ref_strand == "*"
    single <- pairs$ref_mask == "z" & pairs$ref_strand %in% c("+", "-")
    start <- pairs$allele_offset + 1L - (pairs$ref_strand == "-")
    ref_cpg <- substring(toupper(pairs$ref), start, start + 1L) == "CG"
    alt_cpg <- substring(pairs$alt, start, start + 1L) == "CG"
    (dyad | single) & pairs$ref_context == "CG" & start >= 1L & ref_cpg & alt_cpg
}

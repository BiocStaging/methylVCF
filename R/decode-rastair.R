#' Decode Rastair measurements for one sample and record
#'
#' The `CPG` and `CPGnovo` flags determine the reference and de novo slots in
#' that order. Context and the local sequence determine the allele index and
#' strand when the assignment is unique. Measurements are reported at the
#' context level with a missing haplotype.
#'
#' @param record Parsed VCF record with `pos`, `end`, `ref`, `alt` and `info`.
#' @param fields Named list of the selected sample's FORMAT values.
#' @return A call table with one row per reported context on reference bases
#'   and supported biallelic SNVs. Depth discrepancies are recorded by slot
#'   in the `issues` attribute. Inconsistent slot layouts signal an error.
#' @noRd
decode_rastair <- function(record, fields) {
    if (!any(methyl_format_fields %in% names(fields))) { return(empty_calls()) }
    reference_cpg <- isTRUE(record$info[["CPG"]])
    de_novo_cpg <- isTRUE(record$info[["CPGnovo"]])
    origin <- c("reference", "de_novo")[c(reference_cpg, de_novo_cpg)]
    if (!length(origin)) {
        present <- vapply(methyl_format_fields, function(name) {
            field(fields, name) != "."
        }, logical(1L))
        if (any(present)) {
            abort_methyl_record(
                "Rastair methylation values require CPG or CPGnovo flags."
            )
        }
        return(empty_calls())
    }

    alleles <- c(record$ref, record$alt)
    single_bases <- all(alleles %in% c("A", "C", "G", "T"))
    single_site <- length(alleles) <= 2L && record$end == record$pos
    if (!single_bases || !single_site) {
        abort_methyl_record(
            paste(
                "Rastair methylation is supported for reference bases",
                "and biallelic SNVs."
            ),
            unsupported = TRUE
        )
    }

    values <- read_methyl_values(fields)
    combined <- length(origin) == 2L && all(lengths(values) == 1L)
    if (combined && any(!is.na(unlist(values, use.names = FALSE)))) {
        abort_methyl_record(
            paste(
                "Rastair reports a single methylation slot for two CpG",
                "contexts, so the measurements cannot be assigned."
            ),
            unsupported = TRUE
        )
    }
    values <- lapply(
        values, align_methyl_values, n = length(origin),
        reason = paste(
            "Rastair methylation vector length does not match",
            "the reported CpG contexts."
        )
    )
    slots <- rastair_context_alleles(alleles, origin)
    invalid_depth <- which(values$modified_depth > values$informative_depth)
    keep <- setdiff(seq_along(origin), invalid_depth)
    calls <- new_calls(
        position = rep(record$pos, length(keep)),
        strand = slots$strand[keep],
        context = "CG",
        scope = "context",
        fraction = values$fraction[keep],
        informative_depth = values$informative_depth[keep],
        modified_depth = values$modified_depth[keep],
        allele_index = slots$index[keep],
        motif_origin = origin[keep],
        value_index = keep
    )
    if (length(invalid_depth)) {
        attr(calls, "issues") <- data.frame(
            type = "invalid",
            allele_index = slots$index[invalid_depth],
            allele_offset = 0L,
            value_index = invalid_depth,
            reason = "Rastair ADM5mC exceeds DPM5mC."
        )
    }
    calls
}

#' Assign Rastair context slots to locally identifiable alleles
#'
#' @param alleles REF followed by literal ALT bases at the VCF position.
#' @param origin Reported reference and/or de novo CpG slots in caller order.
#' @return Zero-based allele indices and strands. Ambiguous assignments get
#'   a missing index and "*" strand.
#' @noRd
rastair_context_alleles <- function(alleles, origin) {
    strand_by_base <- c(C = "+", G = "-")
    candidates <- which(alleles %in% names(strand_by_base)) - 1L
    reference <- which(origin == "reference")
    if (length(reference) && !0L %in% candidates) {
        abort_methyl_record("Rastair reference CpG requires REF C or G.")
    }
    index <- rep(NA_integer_, length(origin))
    index[reference] <- 0L

    # A reference CpG already assigned here cannot also fill the de novo slot.
    candidates <- setdiff(candidates, index[reference])
    if (length(candidates) == 1L) { index[origin == "de_novo"] <- candidates }
    strand <- unname(strand_by_base[alleles[index + 1L]])
    strand[is.na(strand)] <- "*"
    list(index = index, strand = strand)
}

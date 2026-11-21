#' Decode DRAGEN measurements for one sample and record
#'
#' `INFO/M5mC` masks are indexed in REF/ALT order. FORMAT measurements
#' follow the called genotype slots. A phased homozygote can report a single
#' pooled vector. `Z` marks a CpG measured across both strands. The lowercase
#' masks describe individual cytosines whose strand follows from the allele.
#'
#' @param record Parsed VCF record with `pos`, `end`, `ref`, `alt` and `info`.
#' @param fields Named list of the selected sample's FORMAT values.
#' @return A call table with mappable genomic anchors. Unmappable slots are
#'   recorded in the `issues` attribute with their original allele and
#'   value indices and allele offset. Inconsistent data signal an error.
#' @noRd
decode_dragen <- function(record, fields) {
    if (!any(methyl_format_fields %in% names(fields))) { return(empty_calls()) }
    mask <- field(record$info, "M5mC")
    values <- read_methyl_values(fields)
    missing <- all(is.na(unlist(values, use.names = FALSE)))
    if (identical(mask, ".")) {
        if (missing) { return(empty_calls()) }
        abort_methyl_record( "DRAGEN methylation values require an INFO/M5mC mask")
    }

    genotype <- field(fields, "GT")
    if (missing && grepl("^\\.([/|]\\.)*$", genotype)) { return(empty_calls()) }
    layout <- dragen_allele_layout(record, mask, genotype)
    slots <- dragen_measurement_slots(layout, genotype, values)
    n <- length(slots$offset)
    if (!n && missing) { return(empty_calls()) }

    values <- lapply(
        values, 
        align_methyl_values, 
        n = n,
        reason = "DRAGEN methylation value count does not match INFO/M5mC and GT"
    )

    projection <- project_dragen_slots(record, layout, slots)
    keep <- is.na(projection$reason)
    index <- slots$index[keep]
    contexts <- c(Z = "CG", z = "CG", x = "CHG", h = "CHH")
    calls <- new_calls(
        position = record$pos + slots$offset[keep],
        strand = projection$strand[keep],
        context = unname(contexts[slots$code[keep]]),
        scope = slots$scope,
        fraction = values$fraction[keep],
        informative_depth = values$informative_depth[keep],
        modified_depth = values$modified_depth[keep],
        allele_index = index,
        haplotype = slots$haplotype[keep],
        motif_origin = ifelse(index == 0L, "reference", "alternate"),
        value_index = which(keep),
        allele_offset = slots$offset[keep]
    )

    if (any(!keep)) {
        attr(calls, "issues") <- data.frame(
            type = "unsupported",
            allele_index = slots$index[!keep],
            allele_offset = slots$offset[!keep],
            value_index = which(!keep),
            reason = projection$reason[!keep]
        )
    }
    calls
}

#' Match DRAGEN genotype indices to sequence alleles and validate their masks
#'
#' @param record Parsed VCF record with coordinates and REF/ALT alleles.
#' @param mask The original INFO/M5mC value.
#' @param genotype The original FORMAT/GT value.
#' @return Masks and sequences in VCF order, the called allele indices, and
#'   flags for literal alleles and multibase reference blocks.
#' @noRd
dragen_allele_layout <- function(record, mask, genotype) {
    if (!grepl("^[0-9]+[/|][0-9]+$", genotype)) {
        if (grepl("^([0-9]+|\\.)([/|]([0-9]+|\\.))*$", genotype)) {
            abort_methyl_record(
                "DRAGEN methylation requires a called diploid FORMAT/GT",
                unsupported = TRUE
            )
        }
        abort_methyl_record("invalid DRAGEN FORMAT/GT")
    }

    sequences <- toupper(c(record$ref, record$alt))
    mask_error <- "invalid DRAGEN INFO/M5mC masks for REF and ALT"
    if (!is.character(mask) || !grepl("^[.Zzxh]+(,[.Zzxh]+)*$", mask)) {
        abort_methyl_record(mask_error)
    }
    masks <- strsplit(mask, ",", fixed = TRUE)[[1L]]
    if (length(masks) != length(sequences)) { abort_methyl_record(mask_error) }
    indices <- as.double(strsplit(genotype, "[/|]")[[1L]])
    if (any(indices >= length(sequences))) {
        abort_methyl_record(
            "DRAGEN FORMAT/GT allele is outside the REF and ALT sequences"
        )
    }
    alleles <- as.integer(indices)

    span <- record$end - record$pos + 1L
    block <- span > 1L && nchar(record$ref) == 1L &&
        all(record$alt %in% c("<NON_REF>", "<*>"))
    if (block && any(alleles != 0L)) {
        abort_methyl_record(
            "DRAGEN multi-base blocks require a reference genotype",
            unsupported = TRUE
        )
    }
    literal <- grepl("^[ACGTN]+$", sequences)
    widths <- nchar(sequences)
    if (block) { widths[[1L]] <- span }
    measured <- literal & grepl("[Zzxh]", masks)
    if (any(measured & nchar(masks) != widths)) {
        abort_methyl_record(
            "DRAGEN methylation mask width does not match the allele or block"
        )
    }
    list(
        sequences = sequences, 
        masks = masks, 
        alleles = alleles,
        literal = literal, 
        block = block
    )
}

#' Expand the called DRAGEN masks in the order of the FORMAT measurements
#'
#' @param layout The validated output of dragen_allele_layout().
#' @param genotype Original diploid GT, including its phase separator.
#' @param values Parsed methylation and depth vectors for this sample.
#' @return Allele indices, offsets, mask codes and local haplotype per slot,
#'   together with the scope of the measurements.
#' @noRd
dragen_measurement_slots <- function(layout, genotype, values) {
    alleles <- layout$alleles
    phased <- grepl("|", genotype, fixed = TRUE)
    if (phased && alleles[[1L]] == alleles[[2L]]) {
        # A phased homozygote may report a single pooled vector.
        mask <- layout$masks[[alleles[[1L]] + 1L]]
        sites <- nchar(gsub(".", "", mask, fixed = TRUE))
        phased <- sites > 0L && any(lengths(values) == 2L * sites)
    }
    slots <- if (phased) seq_along(alleles) else which(!duplicated(alleles))
    masks <- strsplit(layout$masks[alleles[slots] + 1L], "", fixed = TRUE)
    offsets <- lapply(masks, function(mask) which(mask != ".") - 1L)
    codes <- lapply(masks, function(mask) mask[mask != "."])
    count <- lengths(offsets)
    scope <- if (phased) {
        "haplotype"
    } else if (length(slots) == 1L) {
        "pooled"
    } else {
        "allele"
    }
    copies <- if (phased) slots else rep(NA_integer_, length(slots))
    list(
        index = rep(alleles[slots], count),
        offset = unlist(offsets, use.names = FALSE),
        code = unlist(codes, use.names = FALSE),
        haplotype = rep(copies, count),
        scope = scope
    )
}

#' Determine the strand and genomic mappability of each DRAGEN slot
#'
#' @param record Parsed VCF record providing the REF sequence.
#' @param layout Validated sequences and reference block information.
#' @param slots Allele indices, offsets and mask codes in measurement order.
#' @return Strand and an omission reason per slot. A missing reason identifies
#'   a mappable slot. Inconsistent strand evidence signals a record error.
#' @noRd
project_dragen_slots <- function(record, layout, slots) {
    index <- slots$index + 1L
    sequence <- layout$sequences[index]
    offset <- slots$offset
    base <- substring(sequence, offset + 1L, offset + 1L)
    single <- slots$code != "Z"

    # Length-changing alleles only map directly at the shared VCF anchor.
    indel <- nchar(sequence) != nchar(record$ref)
    ref_anchor <- substr(layout$sequences[[1L]], 1L, 1L)
    shared_anchor <- offset == 0L & base == ref_anchor
    reason <- rep(NA_character_, length(offset))
    reason[indel & !shared_anchor] <- "DRAGEN indel offset has no unambiguous reference coordinate"
    reason[!layout$literal[index]] <- "DRAGEN measured allele is not a literal sequence"
    reason[single & layout$block & offset > 0L] <- paste(
        "DRAGEN single-cytosine strand inside a block",
        "requires reference sequence"
    )
    if (any(is.na(reason) & single & !base %in% c("C", "G"))) {
        abort_methyl_record(
            "DRAGEN single-cytosine mask does not match a C or G allele"
        )
    }
    strand <- rep("*", length(offset))
    strand[single & base == "C"] <- "+"
    strand[single & base == "G"] <- "-"
    list(strand = strand, reason = reason)
}

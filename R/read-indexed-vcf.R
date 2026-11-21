#' Read indexed VCF data in source order
#'
#' @param file Path to a BGZF-compressed VCF or gVCF.
#' @param index Path to its TBI or CSI index.
#' @param by_contig Reduced query intervals keyed by contig.
#' @param decode_batch Function that decodes lines against the full query.
#' @return A list of decoded call tables in source order.
#' @noRd
read_indexed_batches <- function(file, index, by_contig, decode_batch) {
    fail <- function(error) {
        cli::cli_abort(
            "Could not query {.path {file}} using {.path {index}}.",
            class = "methylVCF_index_error", parent = error, call = NULL
        )
    }
    tabix <- tryCatch(open_indexed_vcf(file, index), error = fail)
    on.exit(close(tabix), add = TRUE)
    contigs <- tryCatch(Rsamtools::seqnamesTabix(tabix), error = fail)

    batches <- list()
    for (chrom in contigs[contigs %in% names(by_contig)]) {
        regions <- by_contig[[chrom]]
        previous_end <- c(0L, regions$end)
        for (i in seq_len(nrow(regions))) {
            range <- GenomicRanges::GRanges(
                chrom, IRanges::IRanges(regions$start[[i]], regions$end[[i]])
            )
            lines <- tryCatch(
                Rsamtools::scanTabix(tabix, param = range)[[1L]],
                error = fail
            )
            if (!length(lines)) { next }
            for (first in seq.int(1L, length(lines), by = 4096L)) {
                last <- min(length(lines), first + 4095L)
                batches[[length(batches) + 1L]] <- decode_batch(
                    lines[seq.int(first, last)], previous_end = previous_end[[i]]
                )
            }
        }
    }
    batches
}

#' Find the VCF index that the regional reader would select
#' @param file Path to the VCF or gVCF.
#' @return Path to the preferred CSI or TBI index, or `NULL` if absent.
#' @noRd
find_vcf_index <- function(file) {
    candidates <- paste0(file, c(".csi", ".tbi"))
    indices <- candidates[file.exists(candidates)]
    if (length(indices)) indices[[1L]]
}

#' Open a BGZF VCF and its index through Rsamtools
#' @param file Path to the VCF or gVCF.
#' @param index Path to the associated CSI or TBI index.
#' @return An open TabixFile. The caller must close it.
#' @noRd
open_indexed_vcf <- function(file, index) {
    tabix <- Rsamtools::TabixFile(file, index = index)
    open(tabix)
    tabix
}

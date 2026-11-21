library(methylVCF)

region <- function(contig, start, end = start) {
    GenomicRanges::GRanges(contig, IRanges::IRanges(start, end))
}

dragen_file <- system.file(
    "extdata", "dragen-example.g.vcf", package = "methylVCF"
)
rastair_file <- system.file(
    "extdata", "rastair-example.vcf", package = "methylVCF"
)
dragen <- readMethylVcf(
    dragen_file, region("chr21", 13204474, 13204534), format = "dragen"
)
rastair <- readMethylVcf(
    rastair_file, region("chr19", 6105750), format = "rastair"
)

columns <- c(
    "seqnames", "start", "strand", "scope", "allele", "haplotype",
    "fraction", "informative_depth", "modified_depth"
)
print(as.data.frame(dragen)[columns], row.names = FALSE)
print(as.data.frame(rastair)[columns], row.names = FALSE)

# Keep the genomic object available, and use a tidy view for reporting.
has_tidy_tools <- requireNamespace("dplyr", quietly = TRUE) &&
    requireNamespace("tibble", quietly = TRUE)
if (has_tidy_tools) {
    assigned <- rastair |>
        as.data.frame() |>
        tibble::as_tibble() |>
        dplyr::filter(!is.na(.data$allele), .data$informative_depth > 0L) |>
        dplyr::select(
            "sample", "seqnames", "start", "allele", "fraction",
            "informative_depth"
        )
    print(assigned)
}

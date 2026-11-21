vcf_region <- function(chrom, start, end = start) {
    GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
}

dragen_file <- system.file(
    "extdata", "dragen-example.g.vcf", package = "methylVCF"
)
rastair_file <- system.file(
    "extdata", "rastair-example.vcf", package = "methylVCF"
)

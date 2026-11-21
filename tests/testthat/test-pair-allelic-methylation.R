read_allelic_fixture <- function(records) {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    writeLines(c(lines[startsWith(lines, "#")], records), file)
    region <- vcf_region("chr21", 100L, 200L)
    readMethylVcf(file, region)
}

test_that("REF/ALT slots follow allele indices and respect strand", {
    region <- vcf_region("chr21", 13204474L, 13204550L)
    pairs <- pairAllelicMethylation(readMethylVcf(dragen_file, region))
    expect_s4_class(pairs, "GRanges")
    expect_identical(GenomicRanges::start(pairs), c(13204520L, 13204534L))
    expect_identical(pairs$pair_status, rep("different_geometry", 2L))
    expect_identical(pairs$genotype, c("1|0", "0|1"))
    expect_identical(pairs$ref_haplotype, c(2L, 1L))
    expect_identical(pairs$alt_haplotype, c(1L, 2L))
    expect_identical(pairs$ref_strand, c("+", "+"))
    expect_identical(pairs$alt_strand, c("-", "-"))
    expect_identical(pairs$ref_fraction, c(0, NA_real_))
    expect_identical(pairs$alt_fraction, c(0.8, 0.7))
    expect_identical(pairs$delta, rep(NA_real_, 2L))
})

test_that("a CpG delta requires a shared motif, geometry and coverage", {
    format <- "GT:M5mC:DPM5mC:PS"
    calls <- read_allelic_fixture(c(
        paste(
            "chr21", 100, ".", "ACGTA", "TCGTG,<NON_REF>", 60, "PASS",
            "M5mC=.Z...,.Z...,.", format, "1|0:0.8,0:9,14:85", sep = "\t"
        ),
        paste(
            "chr21", 110, ".", "ACGTA", "TCGTG,<NON_REF>", 60, "PASS",
            "M5mC=.Z...,.Z...,.", format, "0|1:0,0.8:0,9:85", sep = "\t"
        ),
        paste(
            "chr21", 120, ".", "ACGTA", "TCGTG,<NON_REF>", 60, "PASS",
            "M5mC=.Z...,.z...,.", format, "0|1:0.1,0.8:14,9:85", sep = "\t"
        ),
        paste(
            "chr21", 130, ".", "C", "G,<NON_REF>", 60, "PASS",
            "M5mC=Z,Z,.", format, "0|1:0,0.8:14,9:85", sep = "\t"
        )
    ))
    pairs <- pairAllelicMethylation(calls)
    expect_identical(GenomicRanges::start(pairs), c(101L, 111L, 121L, 130L))
    expect_identical(
        pairs$pair_status, c(
            "comparable", "insufficient_coverage", "different_geometry",
            "unverified_cpg"
        )
    )
    expect_identical(pairs$ref_informative_depth, c(14L, 0L, 14L, 14L))
    expect_identical(pairs$alt_informative_depth, rep(9L, 4L))
    expect_identical(pairs$delta, c(0.8, NA_real_, NA_real_, NA_real_))

    slots <- calls[calls$record_start == 100L]
    expect_identical(
        pairAllelicMethylation(slots[1L])$pair_status,
        "missing_call"
    )
    expect_identical(
        pairAllelicMethylation(slots[c(1L, 1L, 2L)])$pair_status,
        "ambiguous"
    )
    expect_error(
        pairAllelicMethylation(IRanges::shift(slots, 1L)),
        "original one-base anchors"
    )
    slots$informative_depth[1L] <- 0L
    expect_identical(
        pairAllelicMethylation(slots)$pair_status,
        "insufficient_coverage"
    )
})

test_that("single-sided masks, empty output and alternate index are explicit", {
    format <- "GT:M5mC:DPM5mC"
    calls <- read_allelic_fixture(c(
        paste(
            "chr21", 100, ".", "ACGTA", "TCGTG,<NON_REF>", 60, "PASS",
            "M5mC=.Z...,.,.", format, "0/1:0:0", sep = "\t"
        ),
        paste(
            "chr21", 110, ".", "T", "<NON_REF>,C", 60, ".",
            "M5mC=.,.,Z", format, "0/2:0.9:10", sep = "\t"
        ),
        paste(
            "chr21", 120, ".", "T", "C,<NON_REF>", 20, "LowGQ",
            "M5mC=.,Z,.", format, "1/0:0.9:10", sep = "\t"
        ),
        paste(
            "chr21", 130, ".", "T", "C,*", 60, "PASS",
            "M5mC=.,Z,.", format, "0/1:0.9:10", sep = "\t"
        )
    ))
    pairs <- pairAllelicMethylation(calls)
    expect_identical(pairs$pair_status, c("ref_only", "alt_only"))
    expect_identical(pairs$alt_index, c(1L, 2L))
    expect_identical(pairs$alt, c("TCGTG", "C"))
    expect_identical(pairs$alt_mask, c(".", "Z"))
    expect_identical(pairs$ref_mask, c("Z", "."))
    expect_identical(pairs$ref_informative_depth, c(0L, NA_integer_))
    expect_identical(
        pairs$ref_coverage_status,
        c("no_informative_depth", NA_character_)
    )
    expect_identical(pairs$delta, rep(NA_real_, 2L))

    all <- pairAllelicMethylation(calls, filters = NULL)
    expect_length(all, 3L)
    expect_identical(
        S4Vectors::metadata(all)["pair_filters"],
        list(pair_filters = NULL)
    )
    empty <- pairAllelicMethylation(calls[FALSE])
    expect_length(empty, 0L)
    expect_identical(
        names(S4Vectors::mcols(empty)),
        names(S4Vectors::mcols(pairs))
    )
    expect_type(empty$delta, "double")
    expect_true("methylation_mask" %in% names(S4Vectors::mcols(calls)))
    expect_false("pairAllelicMethylation" %in% getNamespaceExports("methylVCF"))
})

test_that("complementary masks from distinct records never form a pair", {
    format <- "GT:M5mC:DPM5mC"
    calls <- read_allelic_fixture(c(
        paste(
            "chr21", 100, ".", "C", "G,<NON_REF>", 60, "PASS",
            "M5mC=Z,.,.", format, "0/1:0.4:10", sep = "\t"
        ),
        paste(
            "chr21", 100, ".", "C", "G,<NON_REF>", 60, "PASS",
            "M5mC=.,Z,.", format, "0/1:0.9:10", sep = "\t"
        )
    ))
    pairs <- pairAllelicMethylation(calls)
    expect_identical(pairs$pair_status, c("ref_only", "alt_only"))
    expect_identical(pairs$delta, rep(NA_real_, 2L))
})

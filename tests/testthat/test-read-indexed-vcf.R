indexed_fixture <- function(file) {
    target <- tempfile("methyl VCF ", fileext = ".vcf.gz")
    Rsamtools::bgzip(file, dest = target)
    Rsamtools::indexTabix(target, format = "vcf")
    target
}

expect_indexed_parity <- function(calls, sequential) {
    metadata <- S4Vectors::metadata(calls)
    expect_identical(metadata$access, "tabix")
    metadata[c("access", "file")] <- S4Vectors::metadata(sequential)[
        c("access", "file")
    ]
    S4Vectors::metadata(calls) <- metadata
    expect_identical(calls, sequential)
}

test_that("TBI and CSI return the same calls as sequential reading", {
    dragen_tbi <- indexed_fixture(dragen_file)
    rastair_tbi <- indexed_fixture(rastair_file)
    temporary <- c(dragen_tbi, rastair_tbi)
    on.exit(unlink(c(temporary, paste0(temporary, ".tbi"))), add = TRUE)
    dragen_csi <- test_path("fixtures", "dragen-example.g.vcf.gz")
    rastair_csi <- test_path("fixtures", "rastair-example.vcf.gz")

    dragen_query <- vcf_region(
        c("chr21", "chrAbsent", "chr21", "chr21", "chr21"),
        c(13204539L, 1L, 13204473L, 13204474L, 13204477L),
        c(13204550L, 1L, 13204474L, 13204474L, 13204477L)
    )
    rastair_query <- vcf_region("chr19", c(6105744L, 6105750L))
    dragen <- readMethylVcf(dragen_file, dragen_query)
    rastair <- readMethylVcf(rastair_file, rastair_query, format = "rastair")
    for (file in c(dragen_tbi, dragen_csi)) {
        expect_indexed_parity(readMethylVcf(file, dragen_query), dragen)
    }
    for (file in c(rastair_tbi, rastair_csi)) {
        expect_indexed_parity(readMethylVcf(file, rastair_query, format = "rastair"), rastair)
    }
    empty <- readMethylVcf(dragen_csi, vcf_region("chrAbsent", 1L))
    expect_length(empty, 0L)
    expect_identical(S4Vectors::metadata(empty)$access, "tabix")
})

test_that("disjoint windows preserve source duplicates and later anchors", {
    lines <- readLines(dragen_file)
    header <- startsWith(lines, "#")
    block <- lines[!header][[1L]]
    plain <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(plain), add = TRUE)
    writeLines(c(lines[header], block, block), plain)
    indexed <- indexed_fixture(plain)
    on.exit(unlink(c(indexed, paste0(indexed, ".tbi"))), add = TRUE)
    query <- vcf_region("chr21", c(13204472L, 13204474L, 13204477L))

    calls <- readMethylVcf(indexed, query)
    expect_indexed_parity(calls, readMethylVcf(plain, query))
    expect_identical(GenomicRanges::start(calls), rep(13204474L, 2L))
})

test_that("an unindexed BGZF falls back and a corrupt index raises an error", {
    file <- tempfile(fileext = ".vcf.gz")
    index <- paste0(file, ".tbi")
    on.exit(unlink(c(file, index)), add = TRUE)
    Rsamtools::bgzip(dragen_file, dest = file)
    query <- vcf_region("chr21", 13204474L)
    calls <- readMethylVcf(file, query)
    expect_length(calls, 1L)
    expect_identical(S4Vectors::metadata(calls)$access, "sequential")

    writeLines("Invalid index", index)
    expect_error(
        readMethylVcf(file, query, on_problem = "report"),
        class = "methylVCF_index_error"
    )
})

test_that("reference blocks preserve source coordinates and metadata", {
    block_region <- vcf_region("chr21", 13204474L)
    block <- readMethylVcf(dragen_file, block_region, format = "dragen")
    block_table <- as.data.frame(block)

    expect_identical(as.character(block_table$seqnames), "chr21")
    expect_identical(as.character(block_table$strand), "*")
    expect_identical(block_table$record_end, 13204512L)
    expect_identical(block_table$modified_depth, NA_integer_)
    expect_identical(block_table$allele_offset, 3L)
    expect_identical(block_table$coverage_status, "observed")
    expect_identical(block_table$sample, "dragen_demo")
    expect_identical(block_table$genotype, "0/0")
    expect_identical(block_table$haplotype, NA_integer_)
    expect_equal(
        Seqinfo::seqlengths(block)[["chr21"]],
        46709983L
    )
    expect_equal(nrow(S4Vectors::metadata(block)$diagnostics), 0L)

    selected <- readMethylVcf(
        dragen_file,
        block_region,
        format = "dragen",
        samples = "dragen_demo"
    )
    expect_identical(as.data.frame(selected), block_table)

    named <- readMethylVcf(
        dragen_file,
        block_region,
        samples = c(alias = "dragen_demo")
    )
    expect_identical(names(named), "alias")
    expect_identical(named$sample, "dragen_demo")
})

test_that("overlapping windows select the reported anchors once", {
    query <- vcf_region("chr21", 13204474L)
    block_table <- as.data.frame(
        readMethylVcf(dragen_file, query, format = "dragen")
    )
    overlapping <- vcf_region(
        "chr21",
        c(13204474L, 13204473L),
        c(13204474L, 13204477L)
    )
    empty <- vcf_region("chr21", 13204475L)

    expect_identical(
        as.data.frame(readMethylVcf(dragen_file, overlapping)),
        block_table
    )
    expect_length(readMethylVcf(dragen_file, empty), 0L)
})

test_that("an empty query returns an empty GRanges and no diagnostics", {
    calls <- readMethylVcf(dragen_file, GenomicRanges::GRanges())
    expect_s4_class(calls, "GRanges")
    expect_length(calls, 0L)
    expect_identical(nrow(S4Vectors::metadata(calls)$diagnostics), 0L)
})

test_that("regional imports retain the declared sequence dictionary", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    contig <- which(startsWith(lines, "##contig="))
    lines[contig] <- '##contig=<ID=chr21,length=46709983,assembly="GRCh38">'
    lines <- append(lines, "##contig=<ID=chrUn,assembly=GRCh38>", after = contig)
    writeLines(lines, file)

    declared <- inspectMethylVcf(file)$seqinfo
    calls <- readMethylVcf(file, vcf_region("chr21", 13204474L))
    empty <- readMethylVcf(file, GenomicRanges::GRanges())
    expect_identical(GenomicRanges::seqinfo(calls), declared)
    expect_identical(GenomicRanges::seqinfo(empty), declared)
    expect_identical(Seqinfo::genome(calls), c(chr21 = "GRCh38", chrUn = "GRCh38"))

    absent <- readMethylVcf(file, vcf_region("chrAbsent", 1L))
    expect_length(absent, 0L)
    expect_identical(Seqinfo::seqlevels(absent), c("chr21", "chrUn", "chrAbsent"))
    expect_identical(Seqinfo::seqlengths(absent)[["chrAbsent"]], NA_integer_)
    expect_identical(Seqinfo::genome(absent)[["chrAbsent"]], NA_character_)
})

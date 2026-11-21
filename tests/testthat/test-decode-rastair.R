test_that("Rastair retains the two strands and their original depth counts", {
    query <- vcf_region("chr19", 6105744L, 6105745L)
    strand_calls <- as.data.frame(
        readMethylVcf(rastair_file, query, format = "rastair")
    )

    expect_identical(strand_calls$start, c(6105744L, 6105745L))
    expect_identical(as.character(strand_calls$strand), c("+", "-"))
    expect_identical(strand_calls$fraction, c(0.8, 0.5))
    expect_identical(strand_calls$informative_depth, c(10L, 14L))
    expect_identical(strand_calls$modified_depth, c(8L, 7L))
})

test_that("Rastair keeps reference and de novo CpG contexts separate", {
    query <- vcf_region("chr19", 6105750L)
    cpg <- as.data.frame(readMethylVcf(rastair_file, query, format = "rastair"))

    expect_identical(cpg$start, rep(6105750L, 2L))
    expect_identical(cpg$genotype, rep("1/1", 2L))
    expect_identical(cpg$haplotype, rep(NA_integer_, 2L))
    expect_identical(cpg$motif_origin, c("reference", "de_novo"))
    expect_identical(cpg$value_index, c(1L, 2L))
    expect_identical(cpg$fraction, c(NA_real_, 1))
    expect_identical(cpg$informative_depth, c(0L, 5L))
    expect_identical(cpg$modified_depth, c(0L, 5L))
    expect_identical(cpg$allele_index[[2L]], 1L)
    expect_identical(cpg$allele[[2L]], "G")
    expect_identical(as.character(cpg$strand[[2L]]), "-")
})

test_that("Rastair reports invalid read counts for each affected context", {
    file <- tempfile(fileext = ".vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(rastair_file)
    records <- paste(
        "chr19", 100:103, ".", "C", c("G", ".", ".", "."), ".", "PASS",
        c("CPG;CPGnovo", rep("CPG", 3L)), "GT:M5mC:DPM5mC:ADM5mC",
        c("0/1:0.8,0.667:10,3:12,2", "0/0:0.5:.:3", "0/0:0.5:4:.", "0/0:0:0:1"),
        sep = "\t"
    )
    writeLines(c(lines[startsWith(lines, "#")], records), file)
    region <- vcf_region("chr19", 100L, 103L)

    error <- expect_error(
        readMethylVcf(file, region, format = "rastair"),
        "Rastair ADM5mC exceeds DPM5mC", class = "methylVCF_import_error"
    )
    expect_identical(error$type, "invalid")
    expect_identical(error$value_index, 1L)

    expect_warning(
        calls <- readMethylVcf(file, region, format = "rastair", on_problem = "report"),
        class = "methylVCF_import_warning"
    )
    issues <- S4Vectors::metadata(calls)$diagnostics
    expect_identical(issues$pos, c(100L, 103L))
    expect_identical(issues$type, rep("invalid", 2L))
    expect_identical(issues$value_index, c(1L, 1L))
    expect_identical(GenomicRanges::start(calls), 100:102)
    expect_identical(calls$motif_origin, c("de_novo", "reference", "reference"))
    expect_identical(calls$value_index, c(2L, 1L, 1L))
    expect_identical(calls$fraction, c(0.667, 0.5, 0.5))
})

test_that("a de novo CpG is assigned only to a unique local allele", {
    lines <- readLines(rastair_file)
    header <- lines[startsWith(lines, "#")]
    records <- c(
        "chr19\t100\t.\tC\tG\t.\tPASS\tCPGnovo\tGT:M5mC:DPM5mC\t0/1:0.6:5",
        "chr19\t101\t.\tC\tT\t.\tPASS\tCPGnovo\tGT:M5mC:DPM5mC\t0/1:0.6:5",
        "chr19\t102\t.\tA\tG\t.\tPASS\tCPGnovo\tGT:M5mC:DPM5mC\t0/1:0.6:5"
    )
    file <- tempfile(fileext = ".vcf")
    on.exit(unlink(file), add = TRUE)
    writeLines(c(header, records), file)
    calls <- readMethylVcf(
        file, vcf_region("chr19", 100L, 102L), format = "rastair"
    )
    expect_identical(calls$allele_index, c(NA_integer_, 0L, 1L))
    expect_identical(
        as.character(GenomicRanges::strand(calls)),
        c("*", "+", "-")
    )
    expect_identical(calls$haplotype, rep(NA_integer_, 3L))
})

test_that("Rastair 2.2.0 retains single slots and reports combined contexts", {
    vcf <- system.file(
        "extdata", "rastair-2.2.0-synthetic.vcf", package = "methylVCF"
    )
    region <- vcf_region("synthetic", 101L, 303L)
    expect_warning(
        calls <- readMethylVcf(
            vcf, region, format = "rastair", on_problem = "report"
        ),
        class = "methylVCF_import_warning"
    )

    metadata <- S4Vectors::metadata(calls)
    expect_identical(metadata$caller_version, "2.2.0")
    expect_identical(metadata$vcf_version, "VCFv4.2")
    expect_identical(GenomicRanges::start(calls), c(101L, 102L, 301L, 303L))
    expect_identical(calls$record_start, c(101L, 102L, 301L, 303L))
    expect_identical(calls$allele, c("C", "G", "C", "G"))
    expect_identical(calls$allele_index, rep(0L, 4L))
    expect_identical(as.character(GenomicRanges::strand(calls)),
        c("+", "-", "+", "-"))
    expect_identical(calls$motif_origin,
        c("reference", "reference", "de_novo", "reference"))
    expect_identical(calls$haplotype, rep(NA_integer_, 4L))
    expect_identical(calls$fraction, c(0.666667, 0.666667, 0.5, 0))
    expect_identical(calls$informative_depth, c(12L, 12L, 12L, 0L))
    expect_identical(calls$modified_depth, c(8L, 8L, 9L, 0L))
    expect_identical(calls$coverage_status[[4L]], "no_informative_depth")
    expect_identical(metadata$diagnostics$pos, 302L)
    expect_identical(metadata$diagnostics$type, "unsupported")
    expect_match(metadata$diagnostics$reason, "single methylation slot")

    error <- expect_error(
        readMethylVcf(vcf, region, format = "rastair"),
        class = "methylVCF_import_error"
    )
    expect_identical(error$pos, 302L)
    expect_identical(error$type, "unsupported")
})

test_that("Rastair distinguishes unmeasured records from invalid CpG annotations", {
    file <- tempfile(fileext = ".vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(rastair_file)
    cases <- rbind(
        c("C", ".", ".", "0/0"),
        c("C", ".", ".", "0/0:.:.:."),
        c("C", ".", ".", "0/0:0.5:10:5"),
        c("CG", "C", "CPG", "0/1:0.5:10:5"),
        c("C", "A,G", "CPG", "0/2:0.5:10:5"),
        c("A", "G", "CPG", "0/1:0.5:10:5"),
        c("C", "G", "CPG;CPGnovo", "0/1:0.7,0.2:10:7,2"),
        c("C", "G", "CPG;CPGnovo", "0/1:.:.:."),
        c("G", ".", "CPG", "0/0:0.5:8:4")
    )
    colnames(cases) <- c("ref", "alt", "info", "sample")
    formats <- c("GT", rep("GT:M5mC:DPM5mC:ADM5mC", nrow(cases) - 1L))
    records <- paste(
        "chr19", 100:108, ".", cases[, "ref"], cases[, "alt"], ".", "PASS",
        cases[, "info"], formats, cases[, "sample"], sep = "\t"
    )
    writeLines(c(lines[startsWith(lines, "#")], records), file)

    expect_warning(
        calls <- readMethylVcf(
            file, vcf_region("chr19", 100L, 108L), format = "rastair", on_problem = "report"
        ),
        class = "methylVCF_import_warning"
    )
    issues <- S4Vectors::metadata(calls)$diagnostics
    expect_identical(issues$pos, 102:106)
    expect_identical(issues$type, c("invalid", "unsupported", "unsupported", "invalid", "invalid"))
    expect_identical(issues$reason, c(
        "Rastair methylation values require CPG or CPGnovo flags.",
        rep("Rastair methylation is supported for reference bases and biallelic SNVs.", 2L),
        "Rastair reference CpG requires REF C or G.",
        "Rastair methylation vector length does not match the reported CpG contexts."
    ))
    expect_identical(GenomicRanges::start(calls), c(107L, 107L, 108L))
    expect_identical(calls$motif_origin, c("reference", "de_novo", "reference"))
    expect_identical(calls$allele_index, c(0L, 1L, 0L))
    expect_equal(calls$fraction, c(NA_real_, NA_real_, 0.5))
    expect_identical(calls$informative_depth, c(NA_integer_, NA_integer_, 8L))
})

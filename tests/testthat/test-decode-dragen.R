write_complex_fixture <- function(fixture) {
    lines <- readLines(dragen_file)
    header <- lines[startsWith(lines, "#")]
    last <- length(header)
    header[last] <- paste(header[last], "second", sep = "\t")
    entries <- rbind(
        c(
            13204600L, "C", "CT,CCA,<NON_REF>", "Z,..,...,.", "0/0:0:0:.",
            "0/0:0.5:8:."
        ),
        c(
            13204610L, "C", "A,G,<NON_REF>", "x,.,x,.",
            "2|0:0,0.6:12,10:13204610", "0/0:0.5:8:."
        ),
        c(
            13204620L, "C", "CC,<NON_REF>", "Z,.Z,.", "1/0:0.9,0.2:10,5:.",
            "0/0:0.3:7:."
        ),
        c(
            13204630L, "C", "G,<NON_REF>", "x,x,.", "0|1:0.2:5,10:13204630",
            "1|0:0.8,0.1:9,11:13204630"
        ),
        c(13204640L, "C", "<NON_REF>", "Z,.", "0/0:.:.:.", "0/0:0.4:.:.")
    )
    records <- cbind(
        "chr21",
        entries[, 1L],
        ".",
        entries[, 2:3],
        ".",
        "PASS",
        paste0("M5mC=", entries[, 4L]),
        "GT:M5mC:DPM5mC:PS",
        entries[, 5:6]
    )
    writeLines(c(header, apply(records, 1L, paste, collapse = "\t")), fixture)
}

test_that("phased methylation slots follow the genotype order", {
    query <- vcf_region("chr21", 13204520L, 13204534L)
    variants <- as.data.frame(
        readMethylVcf(dragen_file, query, format = "dragen")
    )
    alt_first <- variants[variants$start == 13204520L, ]
    ref_first <- variants[variants$start == 13204534L, ]

    expect_identical(alt_first$allele_index, c(1L, 0L))
    expect_identical(as.character(alt_first$allele), c("G", "C"))
    expect_identical(as.character(alt_first$strand), c("-", "+"))
    expect_identical(alt_first$haplotype, c(1L, 2L))
    expect_identical(alt_first$fraction, c(0.8, 0))
    expect_identical(alt_first$informative_depth, c(9L, 14L))
    expect_identical(alt_first$value_index, c(1L, 2L))
    expect_identical(ref_first$allele_index, c(0L, 1L))
    expect_identical(ref_first$fraction, c(NA_real_, 0.7))
    expect_identical(ref_first$informative_depth, c(0L, 10L))
})

test_that("phased homozygotes preserve pooled and per-haplotype measurements", {
    query <- vcf_region("chr21", 13204540L, 13204550L)
    homozygous <- as.data.frame(
        readMethylVcf(dragen_file, query, format = "dragen")
    )

    expect_identical(homozygous$start, c(13204540L, 13204550L, 13204550L))
    expect_identical(homozygous$genotype, rep("1|1", 3L))
    expect_identical(homozygous$phase_set, rep("13204520", 3L))
    expect_identical(homozygous$allele_index, rep(1L, 3L))
    expect_identical(homozygous$allele, rep("C", 3L))
    expect_identical(homozygous$scope, c("pooled", "haplotype", "haplotype"))
    expect_identical(homozygous$haplotype, c(NA_integer_, 1L, 2L))
    expect_identical(homozygous$fraction, c(0.4, 0.2, 0.6))
    expect_identical(homozygous$informative_depth, c(15L, 5L, 10L))
})

test_that("assigned alleles and missing informative depth remain distinct", {
    fixture <- tempfile(fileext = ".vcf")
    on.exit(unlink(fixture), add = TRUE)
    write_complex_fixture(fixture)

    query <- vcf_region("chr21", 13204600L, 13204610L)
    decoded <- as.data.frame(
        readMethylVcf(fixture, query, samples = "dragen_demo")
    )
    expect_identical(decoded$start, c(13204600L, 13204610L, 13204610L))
    expect_identical(decoded$allele, c("C", "G", "C"))
    expect_identical(decoded$allele_index, c(0L, 2L, 0L))
    expect_identical(decoded$allele_offset, rep(0L, 3L))
    expect_identical(decoded$haplotype, c(NA_integer_, 1L, 2L))
    expect_identical(decoded$fraction, c(0, 0, 0.6))
    expect_identical(
        decoded$coverage_status,
        c("no_informative_depth", "observed", "observed")
    )

    unavailable <- as.data.frame(
        readMethylVcf(fixture, vcf_region("chr21", 13204640L))
    )
    expect_identical(unavailable$fraction, c(NA_real_, 0.4))
    expect_identical(unavailable$informative_depth, rep(NA_integer_, 2L))
    expect_identical(
        unavailable$coverage_status,
        c("missing_fraction", "unknown_depth")
    )
})

test_that("report mode collects slot and sample problems in one warning", {
    fixture <- tempfile(fileext = ".vcf")
    on.exit(unlink(fixture), add = TRUE)
    write_complex_fixture(fixture)
    query <- vcf_region("chr21", 13204620L, 13204630L)

    invalid <- expect_error(
        readMethylVcf(fixture, query),
        class = "methylVCF_import_error"
    )
    expect_identical(invalid$chrom, "chr21")
    expect_identical(invalid$pos, 13204620L)
    expect_identical(invalid$sample, "dragen_demo")
    expect_identical(invalid$type, "unsupported")
    expect_identical(invalid$allele_index, 1L)
    expect_identical(invalid$allele_offset, 1L)
    expect_identical(invalid$value_index, 1L)

    warnings <- list()
    partial <- withCallingHandlers(
        readMethylVcf(fixture, query, on_problem = "report"),
        warning = function(w) {
            warnings[[length(warnings) + 1L]] <<- w
            invokeRestart("muffleWarning")
        }
    )
    issues <- S4Vectors::metadata(partial)$diagnostics
    retained <- as.data.frame(partial)

    expect_length(warnings, 1L)
    expect_s3_class(warnings[[1L]], "methylVCF_import_warning")
    expect_identical(warnings[[1L]]$n_issues, 2L)
    expect_identical(issues$pos, c(13204620L, 13204630L))
    expect_identical(issues$type, c("unsupported", "invalid"))
    expect_identical(issues$sample, rep("dragen_demo", 2L))
    expect_identical(issues$allele_index, c(1L, NA_integer_))
    expect_identical(issues$allele_offset, c(1L, NA_integer_))
    expect_identical(issues$value_index, c(1L, NA_integer_))
    expect_identical(issues$reason[[1L]], invalid$reason)
    expect_identical(retained$sample, c("dragen_demo", rep("second", 3L)))
    expect_identical(retained$fraction, c(0.2, 0.3, 0.8, 0.1))
    expect_identical(retained$value_index, c(2L, 1L, 1L, 2L))
})

test_that("DRAGEN skips unmeasured records and diagnoses inconsistent annotations", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    cases <- rbind(
        c("T", ".", "0/1:.:."),
        c("T", "M5mC=Z,.", "./.:.:."),
        c("T", "M5mC=.,.", "0/1:.:."),
        c("T", ".", "0/1:0.5:10"),
        c("T", "M5mC=Z,q", "0/1:0.5:10"),
        c("T", "M5mC=Z", "0/1:0.5:10"),
        c("T", "M5mC=ZZ,.", "0/1:0.5:10"),
        c("T", "M5mC=Z,.", "0/2:0.5:10"),
        c("T", "M5mC=Z,.", "0:0.5:10"),
        c("T", "M5mC=Z,.", "0/x:0.5:10"),
        c("T", "M5mC=.,x", "0/1:0.5:10"),
        c("<NON_REF>", "END=114;M5mC=Z...,.", "0/1:0.5:10"),
        c("T", "M5mC=Z,.", "0/1:0.4:8")
    )
    colnames(cases) <- c("alt", "info", "sample")
    records <- paste(
        "chr21", 100:112, ".", "C", cases[, "alt"], ".", "PASS",
        cases[, "info"], "GT:M5mC:DPM5mC", cases[, "sample"], sep = "\t"
    )
    writeLines(c(lines[startsWith(lines, "#")], records), file)

    expect_warning(
        calls <- readMethylVcf(file, vcf_region("chr21", 100L, 114L), on_problem = "report"),
        class = "methylVCF_import_warning"
    )
    issues <- S4Vectors::metadata(calls)$diagnostics
    expect_identical(issues$pos, 103:111)
    expect_identical(issues$pos[issues$type == "unsupported"], c(108L, 111L))
    expect_identical(
        issues$pos[issues$type == "invalid"], c(103:107, 109:110)
    )
    expect_identical(issues$reason, c(
        "DRAGEN methylation values require an INFO/M5mC mask",
        rep("invalid DRAGEN INFO/M5mC masks for REF and ALT", 2L),
        "DRAGEN methylation mask width does not match the allele or block",
        "DRAGEN FORMAT/GT allele is outside the REF and ALT sequences",
        "DRAGEN methylation requires a called diploid FORMAT/GT",
        "invalid DRAGEN FORMAT/GT",
        "DRAGEN single-cytosine mask does not match a C or G allele",
        "DRAGEN multi-base blocks require a reference genotype"
    ))
    expect_identical(GenomicRanges::start(calls), 112L)
    expect_identical(calls$allele_index, 0L)
    expect_equal(calls$fraction, 0.4)
})

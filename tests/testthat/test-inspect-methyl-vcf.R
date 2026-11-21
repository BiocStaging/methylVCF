test_that("inspection reports caller-specific fields and sample metadata", {
    multisample <- system.file(
        "extdata", "dragen-multisample.g.vcf", package = "methylVCF"
    )
    dragen <- inspectMethylVcf(multisample)

    expect_s3_class(dragen, "methylVCF_inspection")
    expect_identical(dragen$format, "dragen")
    expect_identical(dragen$source, "DRAGEN")
    expect_identical(sub("^VCFv", "", dragen$vcf_version), "4.5")
    expect_identical(dragen$samples, c("sample_a", "sample_b"))
    expect_true(is.na(dragen$reference))
    expect_s4_class(dragen$seqinfo, "Seqinfo")
    expect_identical(Seqinfo::seqlengths(dragen$seqinfo), c(chr21 = 46709983L))

    mask <- with(dragen$fields, section == "INFO" & id == "M5mC")
    fraction <- with(dragen$fields, section == "FORMAT" & id == "M5mC")
    expect_identical(dragen$fields$number[mask], "R")
    expect_identical(dragen$fields$number[fraction], "M")
    expect_identical(dragen$fields$type[fraction], "Float")
    expect_false(any(dragen$checks$status == "error"))

    rastair <- inspectMethylVcf(rastair_file)
    expect_identical(rastair$format, "rastair")
    expect_identical(rastair$caller_version, "2.3.0-rc1")
    expect_identical(rastair$samples, "rastair_demo")
    expect_identical(Seqinfo::seqlengths(rastair$seqinfo), c(chr19 = 61431566L))
    flags <- with(rastair$fields, section == "INFO" & type == "Flag")
    expect_true(all(c("CPG", "CPGnovo") %in% rastair$fields$id[flags]))
    expect_false(any(rastair$checks$status == "error"))
})

test_that("inspection prints declared fields and diagnostic severity", {
    testthat::local_reproducible_output(crayon = TRUE)
    rendered <- capture.output(print(inspectMethylVcf(dragen_file)))
    plain <- trimws(cli::ansi_strip(rendered))
    title <- grep("Methylation VCF inspection", plain, fixed = TRUE)
    expect_identical(plain[title + 1L], "")
    field <- grep("INFO/M5mC", rendered, value = TRUE, fixed = TRUE)
    expect_true(any(cli::ansi_has_any(field)))
    expect_true(any(grepl("Diagnostics", plain, fixed = TRUE)))
    expect_true(paste(cli::symbol$tick, "Header checks: passed") %in% plain)

    testthat::local_reproducible_output()
    mismatch <- inspectMethylVcf(dragen_file, format = "rastair")
    issues <- trimws(format(mismatch))
    expect_false(any(cli::ansi_has_any(issues)))
    expect_true(paste(cli::symbol$cross, "Header checks: failed") %in% issues)
    expect_true(any(startsWith(issues, paste(cli::symbol$warning, "INFO/CPG"))))
})

test_that("missing methylation declarations are reported without reading records", {
    testthat::local_reproducible_output()
    file <- tempfile(fileext = ".vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    header <- lines[startsWith(lines, "#")]
    bad_record <- "this is deliberately not a VCF record"
    writeLines(c(header, bad_record), file)
    expect_false(any(inspectMethylVcf(file)$checks$status == "error"))

    header <- header[!startsWith(header, "##FORMAT=<ID=DPM5mC,")]
    writeLines(c(header, bad_record), file)
    rendered <- format(inspectMethylVcf(file))
    expect_true(paste(cli::symbol$warning, "Header checks: passed with warnings") %in% rendered)

    header <- header[!startsWith(header, "##FORMAT=<ID=M5mC,")]
    writeLines(c(header, bad_record), file)
    info <- inspectMethylVcf(file, format = "dragen")
    expect_identical(info$format, "dragen")
    expect_identical(info$checks$check[info$checks$status == "error"], "FORMAT/M5mC")

    writeLines(header[!startsWith(header, "#CHROM")], file)
    expect_error(inspectMethylVcf(file), "Missing #CHROM", fixed = TRUE)
})

test_that("inspection distinguishes absent, detected, readable and unreadable indexes", {
    testthat::local_reproducible_output()
    directory <- tempfile("inspection fixtures ")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    fixture <- test_path("fixtures", "dragen-example.g.vcf.gz")
    file <- file.path(directory, basename(fixture))
    expect_true(file.copy(fixture, file))

    unindexed <- inspectMethylVcf(file)
    expect_identical(unindexed$index$status, "absent")

    Rsamtools::indexTabix(file, format = "vcf")
    tbi <- inspectMethylVcf(file)
    expect_identical(tbi$index$status, "readable")
    expect_identical(tolower(tbi$index$type), "tbi")
    expect_true("chr21" %in% tbi$index$contigs)

    expect_true(file.copy(paste0(fixture, ".csi"), paste0(file, ".csi")))
    csi <- inspectMethylVcf(file)
    expect_identical(csi$index$status, "readable")
    expect_identical(tolower(csi$index$type), "csi")

    writeLines("Invalid index", paste0(file, ".csi"))
    unprobed <- inspectMethylVcf(file, check_index = FALSE)
    expect_identical(unprobed$index$status, "detected")
    expect_identical(tolower(unprobed$index$type), "csi")

    broken <- inspectMethylVcf(file)
    expect_identical(broken$index$status, "unreadable")
    expect_identical(tolower(broken$index$type), "csi")
    expect_match(broken$index$path, "[.]csi$")
    expect_true(is.character(broken$index$message) && nzchar(broken$index$message))

    detected <- trimws(format(unprobed))
    expect_true(any(startsWith(detected, paste(cli::symbol$info, "Access:"))))
    expect_true(paste(cli::symbol$tick, "Access: TBI index readable") %in% trimws(format(tbi)))
    expect_true(paste(cli::symbol$cross, "Access: CSI index unreadable") %in% trimws(format(broken)))
})

test_that("structured header fields are parsed in any order", {
    file <- tempfile("calls {example}", fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    fraction <- startsWith(lines, "##FORMAT=<ID=M5mC,")
    lines[fraction] <- paste0(
        '##FORMAT=<Number=M,Description="Reported, Number=1, fraction",',
        'ID=M5mC,Type=Float>'
    )
    contig <- startsWith(lines, "##contig=<")
    lines[contig] <- '##contig=<ID=chr21,length=46709983,assembly="GRCh38">'
    writeLines(lines, file)

    info <- inspectMethylVcf(file, format = "rastair")
    expect_identical(info$caller, "dragen")
    expect_identical(info$format, "rastair")
    expect_identical(unname(Seqinfo::genome(info$seqinfo)), "GRCh38")
    status <- info$checks$status
    expect_identical(status[info$checks$check == "caller"], "error")
    expect_identical(status[info$checks$check == "FORMAT/M5mC"], "pass")
    expect_length(readMethylVcf(file, vcf_region("chr21", 13204474L)), 1L)

    rendered <- paste(capture.output(print(info)), collapse = "\n")
    expect_match(rendered, basename(file), fixed = TRUE)
})

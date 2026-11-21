test_that("invalid methylation vectors expose source coordinates and sample", {
    malformed <- tempfile(fileext = ".vcf")
    on.exit(unlink(malformed), add = TRUE)
    previous <- options(cli.condition_width = 30L)
    on.exit(options(previous), add = TRUE)

    lines <- readLines(dragen_file)
    writeLines(
        sub("1|0:0.8,0:9,14:", "1|0:0.8:9,14:", lines, fixed = TRUE),
        malformed
    )

    invalid <- expect_error(
        readMethylVcf(
            malformed,
            vcf_region("chr21", 13204520L),
            format = "dragen"
        ),
        class = "methylVCF_import_error"
    )
    expect_identical(invalid$chrom, "chr21")
    expect_identical(invalid$pos, 13204520L)
    expect_identical(invalid$sample, "dragen_demo")
    expect_identical(invalid$type, "invalid")
    expect_identical(invalid$allele_index, NA_integer_)
    expect_identical(invalid$allele_offset, NA_integer_)
    expect_identical(invalid$value_index, NA_integer_)
    expect_identical(
        invalid$reason,
        "DRAGEN methylation value count does not match INFO/M5mC and GT"
    )
})

test_that("numeric FORMAT fields retain missing values and reject overflow", {
    values <- read_methyl_values(list(M5mC = ".,4.2e-1", DPM5mC = "0,12"))
    expect_identical(values$fraction, c(NA_real_, 0.42))
    expect_identical(values$informative_depth, c(0L, 12L))
    expect_identical(values$modified_depth, NA_integer_)

    error <- expect_error(read_methyl_values(list(DPM5mC = "2147483648")))
    expect_identical(error$reason, "Invalid DPM5mC read depth.")
    error <- expect_error(read_methyl_values(list(M5mC = "1+0i")))
    expect_identical(error$reason, "Invalid M5mC value.")
    error <- expect_error(read_methyl_values(list(M5mC = ".,T")))
    expect_identical(error$reason, "Invalid M5mC value.")
})

test_that("VCF positions and optional contig lengths respect integer bounds", {
    expect_identical(as_position("2147483647"), .Machine$integer.max)
    error <- expect_error(as_position("2147483648"))
    expect_identical(error$reason, "Invalid VCF coordinate.")

    header <- c(
        "<assembly=GRCh38,length=248956422,ID=chr1>", "<ID=chr2>",
        "<ID=chr3,length=2147483648>", "<ID=chr4,length=0>"
    )
    seqinfo <- vcf_contig_seqinfo(header)
    expect_identical(
        Seqinfo::seqlengths(seqinfo),
        c(chr1 = 248956422L, chr2 = NA_integer_,
            chr3 = NA_integer_, chr4 = NA_integer_)
    )
    expect_identical(
        Seqinfo::genome(seqinfo),
        c(chr1 = "GRCh38", chr2 = NA_character_,
            chr3 = NA_character_, chr4 = NA_character_)
    )
})

test_that("a maximum-coordinate SNV is valid and REF overflows are reported", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    header <- lines[startsWith(lines, "#") & !startsWith(lines, "##contig=")]
    position <- .Machine$integer.max
    high <- paste(
        "chr21", position, ".", "C", "T", 60, "PASS",
        "M5mC=Z,.", "GT:M5mC:DPM5mC", "0/1:0.5:10", sep = "\t"
    )
    overflow <- sub("\tC\tT\t", "\tCG\tT\t", high)
    writeLines(c(header, high, overflow), file)
    region <- vcf_region("chr21", position)

    expect_warning(
        calls <- readMethylVcf(file, region, on_problem = "report"),
        class = "methylVCF_import_warning"
    )
    expect_identical(GenomicRanges::start(calls), position)
    expect_identical(calls$record_end, position)
    diagnostic <- S4Vectors::metadata(calls)$diagnostics
    expect_identical(
        diagnostic$reason,
        "Record END exceeds the supported coordinate range."
    )
})

test_that("declared END cannot truncate the REF sequence", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    record <- paste(
        "chr21", 100, ".", "CG", "TG", 60, "PASS",
        "END=100;M5mC=Z.,..", "GT:M5mC:DPM5mC", "0/1:0.5:10", sep = "\t"
    )
    writeLines(c(lines[startsWith(lines, "#")], record), file)
    region <- vcf_region("chr21", 100L)
    expect_error(readMethylVcf(file, region), "END precedes the last REF base")
})

test_that("empty VCF sample names fail at the reader boundary", {
    region <- vcf_region("chr21", 13204474L)
    expect_error(
        readMethylVcf(dragen_file, region, samples = ""),
        "nonempty sample names"
    )
    expect_length(readMethylVcf(dragen_file, region, samples = character()), 0L)

    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    header <- startsWith(lines, "#CHROM")
    lines[header] <- paste0(lines[header], "\t")
    writeLines(lines, file)
    expect_error(readMethylVcf(file, region), "Empty sample name in VCF header")
})

test_that("only matching header declarations update import provenance", {
    lines <- readLines(dragen_file)
    columns <- which(startsWith(lines, "#CHROM"))
    declarations <- c(
        "##source=unit=test", "##source",
        '##DRAGENCommandLine=<ID=dragen,Version="SW: 4.4.6">',
        '##DRAGENCommandLine=<ID=dragen,Version="SW: .">',
        '##FORMAT=<ID=M5mCX,Number=1,Type=String,Description="Other">'
    )
    file <- tempfile(fileext = ".vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- append(lines, declarations, after = columns - 1L)
    writeLines(lines, file)
    query <- vcf_region("chr21", 13204474L)
    metadata <- S4Vectors::metadata(readMethylVcf(file, query))
    expect_identical(metadata$source, "unit=test")
    expect_identical(metadata$caller_version, "4.4.6")
    declared <- !startsWith(lines, "##FORMAT=<ID=M5mC,")
    writeLines(lines[declared], file)
    expect_error(readMethylVcf(file, query), "Missing FORMAT/M5mC")
})

test_that("report mode surfaces unexpected record and decoder errors", {
    region <- vcf_region("chr21", 13204474L)
    fail <- function(...) {
        cli::cli_abort("Unexpected implementation error.", class = "test_fault")
    }
    testthat::with_mocked_bindings(
        expect_error(
            readMethylVcf(dragen_file, region, on_problem = "report"),
            class = "test_fault"
        ),
        decode_dragen = fail
    )
    testthat::with_mocked_bindings(
        expect_error(
            readMethylVcf(dragen_file, region, on_problem = "report"),
            class = "test_fault"
        ),
        parse_vcf_record = fail
    )
})

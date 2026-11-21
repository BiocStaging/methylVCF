test_that("experiment sites follow the declared sequence order", {
    lines <- readLines(dragen_file)
    header <- startsWith(lines, "#")
    block <- lines[!header][[1L]]
    other <- strsplit(block, "\t", fixed = TRUE)[[1L]]
    other[1:2] <- c("chrA", "100")
    other[[8L]] <- sub("END=[0-9]+", "END=141", other[[8L]])
    second <- paste(other, collapse = "\t")
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    writeLines(c(lines[header], block, second), file)

    query <- vcf_region(c("chr21", "chrA"), c(13204474L, 103L))
    calls <- readMethylVcf(file, query)
    experiment <- asMethylExperiment(calls)
    sites <- SummarizedExperiment::rowRanges(experiment)
    expect_identical(
        as.character(GenomicRanges::seqnames(sites)), c("chr21", "chrA")
    )
})

test_that("multisample assays preserve order and distinguish missing values", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)

    lines <- readLines(dragen_file)
    header <- startsWith(lines, "#CHROM")
    records <- !startsWith(lines, "#")
    extra_samples <- c(
        "absent", "covered_zero", "no_depth", "unknown_depth", "no_fraction"
    )
    sample_fields <- c(".", "0/0:0:6", "0/0:0:0", "0/0:0.4:.", "0/0:.:4")
    lines[header] <- paste(
        lines[header],
        paste(extra_samples, collapse = "\t"),
        sep = "\t"
    )
    lines[records] <- paste(
        lines[records],
        paste(sample_fields, collapse = "\t"),
        sep = "\t"
    )
    writeLines(lines, file)

    region <- vcf_region("chr21", 13204474L)
    samples <- append(extra_samples, "dragen_demo", after = 2L)
    calls <- readMethylVcf(file, region, samples = samples)
    experiment <- asMethylExperiment(calls)
    assays <- SummarizedExperiment::assays(experiment, withDimnames = FALSE)
    sample_data <- SummarizedExperiment::colData(experiment)
    imports <- S4Vectors::metadata(experiment)$imports

    expect_identical(colnames(experiment), samples)
    expect_identical(sample_data$vcf_sample, samples)
    expect_identical(sample_data$import_index, rep(1L, 6L))
    expect_identical(imports[[1L]], S4Vectors::metadata(calls))
    expect_identical(assays$fraction[1L, ], c(NA, 0, 0.2381, NA, NA, NA))
    expect_identical(
        assays$reported_fraction[1L, ],
        c(NA, 0, 0.2381, 0, 0.4, NA)
    )
    expect_identical(assays$modified_depth[1L, ], rep(NA_integer_, 6L))
    expect_identical(assays$informative_depth[1L, ], c(NA, 6L, 21L, 0L, NA, 4L))
    expect_identical(
        assays$coverage_status[1L, ],
        c(
            NA, "observed", "observed", "no_informative_depth", "unknown_depth",
            "missing_fraction"
        )
    )
})

test_that("DRAGEN assays keep pooled reference sites and empty samples", {
    query <- vcf_region("chr21", 13204474L, 13204550L)
    variant <- vcf_region("chr21", 13204520L, 13204550L)
    calls <- readMethylVcf(dragen_file, query)
    variants <- readMethylVcf(dragen_file, variant)

    experiment <- asMethylExperiment(list(tumour = calls, normal = variants))
    sites <- SummarizedExperiment::rowRanges(experiment)
    sample_data <- SummarizedExperiment::colData(experiment)
    expect_identical(colnames(experiment), c("tumour", "normal"))
    expect_identical(unname(sample_data$vcf_sample), rep("dragen_demo", 2L))
    expect_identical(sample_data$import_index, 1:2)
    expect_identical(GenomicRanges::start(sites), 13204474L)
    expect_identical(as.character(GenomicRanges::strand(sites)), "*")
    expect_equal(
        unname(SummarizedExperiment::assay(experiment, "fraction")[1L, ]),
        c(0.2381, NA_real_)
    )
    expect_identical(dim(asMethylExperiment(variants)), c(0L, 1L))

    empty <- readMethylVcf(dragen_file, GenomicRanges::GRanges())
    empty_experiment <- asMethylExperiment(empty)
    expect_identical(dim(empty_experiment), c(0L, 1L))
    expect_identical(colnames(empty_experiment), "dragen_demo")

    S4Vectors::metadata(calls)$samples <- NULL
    expect_identical(colnames(asMethylExperiment(calls)), "dragen_demo")
})

test_that("Rastair assays align reference contexts and retain counts", {
    query <- vcf_region("chr19", 6105744L, 6105750L)
    calls <- readMethylVcf(rastair_file, query, format = "rastair")
    first <- calls[GenomicRanges::start(calls) != 6105744L]
    second <- calls[GenomicRanges::start(calls) != 6105745L]
    experiment <- asMethylExperiment(list(first = first, second = second))
    sites <- SummarizedExperiment::rowRanges(experiment)
    assays <- SummarizedExperiment::assays(experiment, withDimnames = FALSE)

    expect_identical(
        GenomicRanges::start(sites),
        c(6105744L, 6105745L, 6105750L)
    )
    expect_identical(
        as.character(GenomicRanges::strand(sites)),
        c("+", "-", "+")
    )
    expect_identical(assays$modified_depth[, 1L], c(NA_integer_, 7L, 0L))
    expect_identical(assays$modified_depth[, 2L], c(8L, NA_integer_, 0L))
    expect_identical(assays$informative_depth[, 2L], c(10L, NA_integer_, 0L))
    expect_identical(assays$reported_fraction[, 2L], c(0.8, NA_real_, NA_real_))
    expect_identical(assays$fraction[, 2L], c(0.8, NA_real_, NA_real_))
    expect_identical(
        assays$coverage_status[, 2L],
        c("observed", NA, "no_informative_depth")
    )
})

test_that("the same anchor can have distinct strands and contexts", {
    query <- vcf_region("chr21", 13204474L)
    calls <- readMethylVcf(dragen_file, query)[rep(1L, 3L)]
    GenomicRanges::strand(calls) <- c("*", "+", "+")
    calls$context <- c("CG", "CG", "CHG")
    sites <- SummarizedExperiment::rowRanges(asMethylExperiment(calls))

    expect_length(sites, 3L)
    expect_identical(S4Vectors::mcols(sites)$context, c("CG", "CG", "CHG"))
    expect_identical(
        as.character(GenomicRanges::strand(sites)),
        c("*", "+", "+")
    )
})

test_that("the default FILTER selection excludes failed overlapping records", {
    query <- vcf_region("chr21", 13204474L)
    calls <- readMethylVcf(dragen_file, query)[rep(1L, 2L)]
    calls$filter <- c("LowGQ", "PASS")
    calls$fraction <- c(0.6, 0.2381)
    experiment <- asMethylExperiment(calls)

    expect_identical(S4Vectors::metadata(experiment)$filters, c("PASS", "."))
    expect_equal(
        SummarizedExperiment::assay(experiment, "fraction")[1L, 1L],
        0.2381
    )
    expect_error(
        asMethylExperiment(calls, filters = NULL),
        class = "methylVCF_duplicate_site"
    )
    expect_error(
        asMethylExperiment(calls, filters = NA_character_),
        "FILTER values"
    )

    calls$filter[2L] <- "."
    expect_identical(dim(asMethylExperiment(calls)), c(1L, 1L))
    expect_identical(
        dim(asMethylExperiment(calls, filters = "PASS")),
        c(0L, 1L)
    )
})

test_that("experiments reject duplicate cells and incompatible references", {
    query <- vcf_region("chr21", 13204474L)
    calls <- readMethylVcf(dragen_file, query)
    expect_error(
        asMethylExperiment(calls[c(1L, 1L)]),
        class = "methylVCF_duplicate_site"
    )
    expect_error(
        asMethylExperiment(list(calls)),
        "unique, nonempty sample names"
    )

    other <- calls
    S4Vectors::metadata(calls)$reference <- "reference-A"
    S4Vectors::metadata(other)$reference <- "reference-B"
    expect_error(
        asMethylExperiment(list(a = calls, b = other)),
        "different VCF references"
    )

    S4Vectors::metadata(other)$format <- "rastair"
    expect_error(asMethylExperiment(list(a = calls, b = other)), "one format")

    other <- calls
    coordinates <- data.frame(
        seqnames = "chr21",
        start = 13204474L,
        end = 13204474L
    )
    shorter_contig <- GenomicRanges::makeGRangesFromDataFrame(
        coordinates, seqinfo = c(chr21 = 46709982L)
    )
    GenomicRanges::seqinfo(other) <- GenomicRanges::seqinfo(shorter_contig)
    expect_error(
        asMethylExperiment(list(a = calls, b = other)),
        "incompatible sequence information"
    )
})

test_that("experiments compare declared assemblies across VCF imports", {
    file <- tempfile(fileext = ".g.vcf")
    on.exit(unlink(file), add = TRUE)
    lines <- readLines(dragen_file)
    contig <- startsWith(lines, "##contig=")
    query <- vcf_region("chr21", 13204474L)
    import <- function(assembly, reference = NULL) {
        contents <- lines
        contents[contig] <- paste0(
            "##contig=<ID=chr21,length=46709983,assembly=", assembly, ">"
        )
        if (!is.null(reference)) {
            contents <- append(contents, paste0("##reference=", reference), after = 1L)
        }
        writeLines(contents, file)
        readMethylVcf(file, query)
    }

    expect_error(
        asMethylExperiment(list(a = import("build_A"), b = import("build_B"))),
        "incompatible genome assemblies"
    )

    a <- import("build_A", "file:///first/reference.fa")
    b <- import("build_A", "file:///second/reference.fa")
    experiment <- asMethylExperiment(list(a = a, b = b))
    expect_identical(dim(experiment), c(1L, 2L))
    expect_identical(Seqinfo::genome(experiment), c(chr21 = "build_A"))
    expect_identical(
        S4Vectors::metadata(experiment)$imports[[2L]]$reference,
        "file:///second/reference.fa"
    )

    unknown <- import(".", "file:///second/reference.fa")
    expect_error(asMethylExperiment(list(a = a, b = unknown)), "different VCF references")
})

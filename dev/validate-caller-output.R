# Compare caller fields with the installed package. See dev/caller-validation.md.
library(methylVCF)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
    cli::cli_abort(
        paste(
            "Usage: Rscript --vanilla dev/validate-caller-output.R",
            "<rastair-reheader.vcf> <HG002.gvcf.gz>"
        )
    )
}
rastair_file <- arguments[[1L]]
hg002_file <- arguments[[2L]]
region <- function(chrom, start, end = start) {
    GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
}

format_fields <- function(record) {
    keys <- strsplit(record[[9L]], ":", fixed = TRUE)[[1L]]
    values <- strsplit(record[[10L]], ":", fixed = TRUE)[[1L]]
    stats::setNames(values, keys)
}

raw_records <- function(lines) {
    strsplit(lines[!startsWith(lines, "#")], "\t", fixed = TRUE)
}

raw_rastair <- function(record) {
    stopifnot("CPG" %in% strsplit(record[[8L]], ";", fixed = TRUE)[[1L]])
    fields <- format_fields(record)
    data.frame(
        position = as.integer(record[[2L]]),
        ref = record[[4L]],
        alt = record[[5L]],
        filter = record[[7L]],
        genotype = unname(fields[["GT"]]),
        fraction = as.numeric(fields[["M5mC"]]),
        informative_depth = as.integer(fields[["DPM5mC"]]),
        modified_depth = as.integer(fields[["ADM5mC"]])
    )
}

rastair_region <- region("chr19", 6103156L, 6143156L)
source <- raw_records(readLines(rastair_file))
expected <- do.call(rbind, lapply(source, raw_rastair))
expected$alt[expected$alt == "."] <- ""
calls <- readMethylVcf(rastair_file, rastair_region, format = "rastair")
metadata <- S4Vectors::metadata(calls)
observed <- as.data.frame(calls)
expected_strand <- unname(c(C = "+", G = "-")[expected$ref])
expected_status <- ifelse(
    expected$informative_depth > 0L, "observed", "no_informative_depth"
)
stopifnot(
    nrow(expected) == 1901L,
    metadata$caller_version == "2.2.0",
    nrow(metadata$diagnostics) == 0L,
    all(as.character(observed$seqnames) == "chr19"),
    identical(observed$start, expected$position),
    identical(observed$end, expected$position),
    identical(observed$record_start, expected$position),
    identical(observed$record_end, expected$position),
    identical(observed$ref, expected$ref),
    identical(observed$alt, expected$alt),
    identical(observed$filter, expected$filter),
    identical(observed$genotype, expected$genotype),
    identical(observed$fraction, expected$fraction),
    identical(observed$informative_depth, expected$informative_depth),
    identical(observed$modified_depth, expected$modified_depth),
    identical(observed$coverage_status, expected_status),
    identical(as.character(observed$strand), expected_strand),
    identical(observed$allele, expected$ref),
    all(observed$allele_index == 0L),
    all(observed$motif_origin == "reference"),
    all(is.na(observed$haplotype))
)
kept <- which(expected$filter == "PASS")
stopifnot(!anyDuplicated(expected$position[kept]))
kept <- kept[order(expected$position[kept])]
experiment <- asMethylExperiment(calls)
sites <- SummarizedExperiment::rowRanges(experiment)
assays <- SummarizedExperiment::assays(experiment, withDimnames = FALSE)
sample_info <- SummarizedExperiment::colData(experiment)
raw_fraction <- expected$fraction[kept]
analysis_fraction <- raw_fraction
analysis_fraction[expected$informative_depth[kept] == 0L] <- NA_real_
stopifnot(
    length(sites) == 1444L,
    identical(GenomicRanges::start(sites), expected$position[kept]),
    identical(assays$fraction[, 1L], analysis_fraction),
    identical(assays$reported_fraction[, 1L], raw_fraction),
    identical(assays$informative_depth[, 1L], expected$informative_depth[kept]),
    identical(assays$modified_depth[, 1L], expected$modified_depth[kept]),
    identical(assays$coverage_status[, 1L], expected_status[kept]),
    identical(colnames(experiment), metadata$samples),
    identical(sample_info$vcf_sample, metadata$samples),
    identical(sample_info$import_index, 1L)
)
cli::cli_inform(
    "Rastair 2.2.0 matched {length(calls)} calls and {length(sites)} sites."
)
print(table(strand = expected_strand, coverage = expected_status))

hg_records <- function(positions) {
    con <- gzfile(hg002_file, open = "rt")
    on.exit(close(con), add = TRUE)
    records <- list()
    repeat {
        lines <- readLines(con, n = 4096L)
        if (!length(lines)) { break }
        chr1 <- raw_records(lines[startsWith(lines, "chr1\t")])
        pos <- as.integer(vapply(chr1, `[[`, character(1L), 2L))
        records <- c(records, chr1[pos %in% positions])
        if (any(pos > max(positions))) { break }
    }
    records
}

anchors <- c(10493L, 10816L, 28692L, 94011L, 171535L)
pos <- c(10486L, anchors[-1L])
slots <- c(2L, 1L, 1L, 1L, 1L)
allele <- c(0L, 0L, 0L, 1L, 1L)
phase <- c(NA_integer_, NA_integer_, 1L, NA_integer_, NA_integer_)
scope <- c("pooled", "pooled", "haplotype", "allele", "pooled")
insertion <- 3042652L
records <- hg_records(c(pos, insertion))
source_pos <- as.integer(vapply(records, `[[`, character(1L), 2L))
stopifnot(identical(source_pos, c(pos, insertion)))
hg_calls <- readMethylVcf(hg002_file, region("chr1", anchors))
stopifnot(length(hg_calls) == length(anchors))
for (i in seq_along(anchors)) {
    raw <- records[[i]]
    fields <- format_fields(raw)
    fraction <- as.numeric(strsplit(fields[["M5mC"]], ",")[[1L]])
    depth <- as.integer(strsplit(fields[["DPM5mC"]], ",")[[1L]])
    stopifnot(
        GenomicRanges::start(hg_calls)[[i]] == anchors[[i]],
        hg_calls$record_start[[i]] == pos[[i]],
        hg_calls$allele_index[[i]] == allele[[i]],
        identical(hg_calls$haplotype[[i]], phase[[i]]),
        hg_calls$scope[[i]] == scope[[i]],
        hg_calls$genotype[[i]] == fields[["GT"]],
        hg_calls$fraction[[i]] == fraction[slots[[i]]],
        hg_calls$informative_depth[[i]] == depth[slots[[i]]]
    )
}
block_info <- strsplit(records[[1L]][[8L]], ";", fixed = TRUE)[[1L]]
block_end <- as.integer(sub("END=", "", grep("^END=", block_info, value = TRUE)))
block_mask <- sub("M5mC=", "", grep("^M5mC=", block_info, value = TRUE))
ref_mask <- strsplit(block_mask, ",", fixed = TRUE)[[1L]][[1L]]
active_offsets <- gregexpr("[Zzxh]", ref_mask)[[1L]] - 1L
stopifnot(
    hg_calls$record_end[[1L]] == block_end,
    pos[[1L]] + active_offsets[slots[[1L]]] == anchors[[1L]]
)

raw <- records[[length(records)]]
sample <- format_fields(raw)
fractions <- as.numeric(strsplit(sample[["M5mC"]], ",")[[1L]])
depth <- as.integer(strsplit(sample[["DPM5mC"]], ",")[[1L]])
info <- strsplit(raw[[8L]], ";", fixed = TRUE)[[1L]]
mask <- sub("M5mC=", "", grep("^M5mC=", info, value = TRUE))
alt3 <- strsplit(mask, ",", fixed = TRUE)[[1L]][[4L]]
offsets <- gregexpr("[Zzxh]", alt3)[[1L]] - 1L
insert_calls <- withCallingHandlers(
    readMethylVcf(hg002_file, region("chr1", insertion), on_problem = "report"),
    methylVCF_import_warning = function(cnd) invokeRestart("muffleWarning")
)
omitted <- S4Vectors::metadata(insert_calls)$diagnostics
stopifnot(
    identical(GenomicRanges::start(insert_calls), rep(insertion, 2L)),
    identical(insert_calls$allele_index, c(0L, 3L)),
    identical(insert_calls$genotype, rep(sample[["GT"]], 2L)),
    identical(insert_calls$fraction, fractions[1:2]),
    identical(insert_calls$informative_depth, depth[1:2]),
    identical(omitted$pos, rep(insertion, 3L)),
    identical(omitted$type, rep("unsupported", 3L)),
    identical(omitted$allele_index, rep(3L, 3L)),
    identical(omitted$allele_offset, offsets[-1L]),
    identical(omitted$value_index, 3:5)
)
cli::cli_inform(
    "HG002 matched five anchors and the insertion with offsets 5, 19 and 31."
)

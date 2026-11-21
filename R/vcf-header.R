#' Read and validate the methylation VCF header
#'
#' @param con Open text connection positioned at the start of the VCF.
#' @param require_methylation Whether to require a FORMAT/M5mC declaration.
#' @return A list with samples, caller evidence, reference, sequence information
#'   and raw INFO and FORMAT declarations. Leaves the connection after
#'   `#CHROM`.
#' @noRd
read_vcf_header <- function(con, require_methylation = TRUE) {
    version <- source <- caller_version <- reference <- NA_character_
    contigs <- info <- format <- character()
    caller <- NA_character_
    dragen_version <- '^<ID=dragen,Version="SW: ([0-9][0-9.]*).*'
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (!length(line) || !startsWith(line, "##")) { break }
        key <- sub("=.*$", "=", line)
        value <- substring(line, nchar(key) + 1L)
        switch(key,
            "##fileformat=" = version <- value,
            "##source=" = source <- value,
            "##reference=" = reference <- value,
            "##rastairVersion=" = {
                caller <- "rastair"
                caller_version <- value
            },
            "##DRAGENCommandLine=" = {
                caller <- "dragen"
                if (grepl(dragen_version, value)) {
                    caller_version <- sub(dragen_version, "\\1", value)
                }
            },
            "##contig=" = {
                if (startsWith(value, "<")) { contigs <- c(contigs, value) }
            },
            "##INFO=" = info <- c(info, value),
            "##FORMAT=" = format <- c(format, value)
        )
    }
    if (!length(line) || !startsWith(line, "#CHROM\t")) {
        cli::cli_abort("Missing #CHROM header in VCF.", call = NULL)
    }
    if (is.na(version) || !grepl("^VCFv[0-9]+\\.[0-9]+$", version)) {
        cli::cli_abort("Missing or invalid VCF fileformat header.", call = NULL)
    }
    samples <- vcf_header_samples(line)
    if (require_methylation && !has_methylation_declaration(format)) {
        cli::cli_abort("Missing FORMAT/M5mC declaration in VCF.", call = NULL)
    }
    if (is.na(caller)) {
        source_caller <- tolower(sub(
            '^[[:space:]"]*(dragen|rastair).*', "\\1", source, ignore.case = TRUE
        ))
        if (source_caller %in% c("dragen", "rastair")) { caller <- source_caller }
    }
    list(
        samples = samples, 
        source = source, 
        version = version,
        caller = caller,
        caller_version = caller_version, 
        reference = reference,
        seqinfo = vcf_contig_seqinfo(contigs),
        info = info, format = format
    )
}

#' Read sample IDs from the VCF column header
#'
#' @param line The tab-delimited column header beginning with #CHROM.
#' @return Unique sample IDs in file order. Malformed columns raise an error.
#' @noRd
vcf_header_samples <- function(line) {
    fixed <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    columns <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    has_sample <- length(columns) > length(fixed)
    standard <- identical(columns[seq_along(fixed)], fixed)
    if (!has_sample || !standard) {
        cli::cli_abort(
            "The VCF must have standard columns and at least one sample.",
            call = NULL
        )
    }
    samples <- columns[-seq_along(fixed)]
    if (endsWith(line, "\t") || any(!nzchar(samples))) {
        cli::cli_abort("Empty sample name in VCF header.", call = NULL)
    }
    if (anyDuplicated(samples)) {
        cli::cli_abort("Duplicate sample names in VCF.", call = NULL)
    }
    samples
}

#' Extract declared contigs, lengths and assemblies
#'
#' @param contigs Character vector of structured VCF `##contig` values.
#' @return A `Seqinfo` in header order, using `NA` for missing assemblies and
#'   missing, zero or unrepresentable lengths. Duplicate IDs signal an error.
#' @noRd
vcf_contig_seqinfo <- function(contigs) {
    declarations <- vcf_declaration_columns(contigs, c("ID", "length", "assembly"))
    has_id <- !is.na(declarations$ID) & nzchar(declarations$ID)
    ids <- declarations$ID[has_id]
    if (anyDuplicated(ids)) {
        cli::cli_abort("Duplicate contig names in VCF header.", call = NULL)
    }
    declared_lengths <- declarations$length[has_id]
    declared <- grep("^[0-9]+$", declared_lengths)
    numbers <- as.double(declared_lengths[declared])
    valid <- numbers >= 1 & numbers <= .Machine$integer.max
    lengths <- rep(NA_integer_, length(ids))
    lengths[declared[valid]] <- as.integer(numbers[valid])
    assemblies <- declarations$assembly[has_id]
    assemblies[assemblies %in% c("", ".")] <- NA_character_
    Seqinfo::Seqinfo(ids, seqlengths = lengths, genome = assemblies)
}

#' Check for the methylation fraction declaration required by the reader
#' @param lines Raw FORMAT declarations from the VCF header.
#' @return Whether FORMAT/M5mC is declared.
#' @noRd
has_methylation_declaration <- function(lines) {
    "M5mC" %in% vcf_declaration_columns(lines, "ID")$ID
}

#' Extract named attributes from structured header declarations
#'
#' @param lines Values delimited by angle brackets in the header.
#' @param keys Names of attributes to extract.
#' @return A named list of character vectors aligned with the input lines.
#' @noRd
vcf_declaration_columns <- function(lines, keys) {
    # Match quoted attributes as a unit, including escaped quotes and commas.
    lines <- gsub("^<|>$", "", lines)
    pattern <- r"((^|,)[^,=]+=("[^"\\]*(?:\\.[^"\\]*)*"|[^,]*))"
    tokens <- regmatches(lines, gregexpr(pattern, lines, perl = TRUE))
    attributes <- function(key) {
        prefix <- paste0(key, "=")
        vapply(tokens, function(record) {
            record <- trimws(sub("^,", "", record))
            value <- record[startsWith(record, prefix)]
            if (!length(value)) { return(NA_character_) }
            sub('^"(.*)"$', "\\1", substring(value[[1L]], nchar(prefix) + 1L))
        }, character(1L))
    }
    stats::setNames(lapply(keys, attributes), keys)
}

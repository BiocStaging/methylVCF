#' Inspect a methylation VCF before importing it
#'
#' Read the header of a local DRAGEN or Rastair VCF and optionally check
#' whether its index opens. The inspection stops at the column header and
#' does not decode methylation records.
#'
#' @param file Path to a local VCF or gVCF, optionally compressed.
#' @param format The expected caller, `"dragen"` or `"rastair"`. With `NULL`,
#'   identify the caller from the header when possible.
#' @param check_index Logical. Try to open the index that [readMethylVcf()]
#'   would select. Set to `FALSE` to check only for its presence.
#'
#' @return A list of class `methylVCF_inspection` with a console summary and
#'   the following components:
#'   - `file`, `caller`, `source`, `vcf_version`, `caller_version` and
#'     `reference` describe the input and its declared provenance. Unavailable
#'     values are `NA`. `format` is the caller used for the header checks.
#'   - `samples` contains VCF sample IDs in file order. `seqinfo` holds
#'     the contigs, lengths and any assembly tags declared per contig.
#'     The file-level reference declaration is returned as `reference`.
#'   - `fields` is a data frame of INFO and FORMAT declarations with
#'     `section`, `id`, `number` and `type` columns.
#'   - `checks` is a data frame of header findings with `check`, `status`
#'     (`"pass"`, `"warning"` or `"error"`) and `message` columns.
#'   - `index` is a list with the chosen `path`, `type` (`"TBI"` or `"CSI"`),
#'     `status` (`"absent"`, `"detected"`, `"readable"` or `"unreadable"`),
#'     indexed `contigs`, and a `message` if opening fails.
#'
#' @details The selected CSI index takes precedence over TBI, including
#'   when the CSI index is unreadable. A readable index can be opened and
#'   its contigs listed. Its relationship to the VCF and the VCF records
#'   themselves remain to be checked during a regional import.
#'
#'   A missing methylation declaration is reported in `checks`. A path
#'   that cannot be opened or a structurally invalid VCF header raises an
#'   error. The checks do not infer missing reference builds or require
#'   a particular declared methylation-vector cardinality.
#'
#' @seealso [readMethylVcf()] to query and interpret the records.
#' @examples
#' file <- system.file("extdata", "dragen-multisample.g.vcf", package = "methylVCF")
#' info <- inspectMethylVcf(file)
#' info
#' info$samples
#' @export
#' @md
inspectMethylVcf <- function(file, format = NULL, check_index = TRUE) {
    validate_vcf_file(file)
    if (!is.logical(check_index) || length(check_index) != 1L || is.na(check_index)) {
        cli::cli_abort("{.arg check_index} must be TRUE or FALSE.", call = NULL)
    }

    con <- open_vcf_text(file)
    on.exit(close(con), add = TRUE)
    header <- read_vcf_header(con, require_methylation = FALSE)
    format <- if (is.null(format)) {
        header$caller
    } else {
        match.arg(format, c("dragen", "rastair"))
    }

    fields <- rbind(
        vcf_field_table(header$info, "INFO"),
        vcf_field_table(header$format, "FORMAT")
    )
    structure(
        list(
            file = file,
            format = format,
            caller = header$caller,
            source = header$source,
            vcf_version = header$version,
            caller_version = header$caller_version,
            reference = header$reference,
            samples = header$samples,
            seqinfo = header$seqinfo,
            fields = fields,
            checks = inspect_header_fields(header, format, fields),
            index = inspect_vcf_index(file, check_index)
        ),
        class = "methylVCF_inspection"
    )
}

#' Tabulate the key attributes from INFO or FORMAT declarations
#'
#' @param lines Structured values from `##INFO` or `##FORMAT` header lines.
#' @param section Header section associated with the declarations.
#' @return A data frame of field IDs, declared cardinalities and types.
#' @noRd
vcf_field_table <- function(lines, section) {
    columns <- vcf_declaration_columns(lines, c("ID", "Number", "Type"))
    data.frame(
        section = rep(section, length(lines)),
        id = columns$ID,
        number = columns$Number,
        type = columns$Type
    )
}

#' Summarize caller and methylation declarations
#'
#' @param header Parsed VCF header.
#' @param format Caller used for the checks, or `NA`.
#' @param fields Tabulated field declarations.
#' @return A data frame of check names, statuses and messages.
#' @noRd
inspect_header_fields <- function(header, format, fields) {
    finding <- function(check, status, message) {
        data.frame(check = check, status = status, message = message)
    }
    caller <- if (is.na(format)) {
        finding("caller", "warning", "Caller not identified. Supply format explicitly.")
    } else if (!is.na(header$caller) && header$caller != format) {
        finding("caller", "error", paste("Selected format differs from header caller", header$caller))
    } else {
        finding("caller", "pass", paste("Checking", format, "declarations"))
    }

    fraction <- if (has_methylation_declaration(header$format)) {
        finding("FORMAT/M5mC", "pass", "Methylation fraction is declared.")
    } else {
        finding("FORMAT/M5mC", "error", "Missing FORMAT/M5mC declaration.")
    }

    declared <- paste(fields$section, fields$id, sep = "/")
    depth <- if ("FORMAT/DPM5mC" %in% declared) {
        finding("FORMAT/DPM5mC", "pass", "Informative depth is declared.")
    } else {
        finding("FORMAT/DPM5mC", "warning", "Informative depth is not declared.")
    }

    context_fields <- switch(
        format,
        dragen = "INFO/M5mC",
        rastair = c("INFO/CPG", "INFO/CPGnovo")
    )

    context <- if (length(context_fields)) {
        label <- paste(context_fields, collapse = " or ")
        if (any(context_fields %in% declared)) {
            finding(label, "pass", "Caller-specific context is declared.")
        } else {
            finding(label, "warning", "Caller-specific context is not declared.")
        }
    }

    rbind(caller, fraction, depth, context)
}

#' Inspect the companion index without querying VCF records
#'
#' @param file Path to a local VCF or gVCF.
#' @param check_index Whether to open the index.
#' @return A list describing the selected index and the outcome of the probe.
#' @noRd
inspect_vcf_index <- function(file, check_index) {
    index <- find_vcf_index(file)
    if (is.null(index)) {
        return(
            list(
                path = NA_character_,
                type = NA_character_,
                status = "absent",
                contigs = character(),
                message = NA_character_
            )
        )
    }

    result <- list(
        path = index,
        type = toupper(sub("^.*[.]", "", index)),
        status = "detected",
        contigs = character(),
        message = NA_character_
    )

    if (!check_index) { return(result) }

    probe <- tryCatch({
        tabix <- open_indexed_vcf(file, index)
        on.exit(close(tabix), add = TRUE)
        Rsamtools::seqnamesTabix(tabix)
    }, error = identity)

    if (inherits(probe, "error")) {
        result$status <- "unreadable"
        result$message <- conditionMessage(probe)
    } else {
        result$status <- "readable"
        result$contigs <- probe
    }
    result
}

#' Print a compact methylation VCF inspection
#' @param x Inspection from [inspectMethylVcf()].
#' @param ... Unused.
#' @return Invisibly, the original inspection.
#' @export
#' @noRd
print.methylVCF_inspection <- function(x, ...) {
    cat(format(x, ...), sep = "\n")
    invisible(x)
}

#' Format a compact methylation VCF inspection
#' @param x Inspection from [inspectMethylVcf()].
#' @param ... Unused.
#' @return A character vector with one line per element.
#' @export
#' @noRd
format.methylVCF_inspection <- function(x, ...) {
    cli::cli_format_method({
        cli::cli_div(theme = list(
            h1 = list("margin-bottom" = 1),
            h2 = list("margin-bottom" = 1),
            "span.vcf-id" = list(color = "cyan")
        ))
        cli::cli_h1("Methylation VCF inspection")
        caller <- if (is.na(x$caller)) "not identified" else {
            switch(x$caller, dragen = "DRAGEN", rastair = "Rastair")
        }
        version <- if (is.na(x$caller_version)) "version not declared" else x$caller_version
        reference <- if (is.na(x$reference)) "not declared" else x$reference

        lengths <- Seqinfo::seqlengths(x$seqinfo)
        contigs <- names(lengths)
        known <- !is.na(lengths) & length(contigs) <= 3L
        contigs[known] <- paste0(
            contigs[known], " (", format(lengths[known], big.mark = ",", trim = TRUE), " bp)"
        )
        cli::cli_dl(
            list(
                File = "{.file {basename(x$file)}}",
                Caller = "{.strong {caller}} ({version})",
                "VCF version" = "{.val {sub('^VCFv', '', x$vcf_version)}}",
                Reference = "{reference}",
                "Sample IDs" = "{cli::col_cyan(inspection_preview(x$samples))}",
                Contigs = "{inspection_preview(contigs)}"
            ),
            labels = c(
                "File",
                "Caller",
                "VCF version",
                "Reference",
                paste0("Samples (", length(x$samples), ")"),
                paste0("Contigs (", length(contigs), ")")
            )
        )

        declared <- paste(x$fields$section, x$fields$id, sep = "/")
        methylation <- c(
            "INFO/M5mC" = "Allele masks", "FORMAT/M5mC" = "Methylation fraction",
            "FORMAT/DPM5mC" = "Informative depth", "FORMAT/ADM5mC" = "Modified depth",
            "INFO/CPG" = "Reference CpG context", "INFO/CPGnovo" = "De novo CpG context"
        )
        displayed <- methylation[names(methylation) %in% declared]
        if (length(displayed)) {
            cli::cli_h2("Declared methylation fields")
            cli::cli_dl(
                as.list(displayed),
                labels = paste0("{.vcf-id ", names(displayed), "}"),
                class = "inspection-indent"
            )
        }

        cli::cli_h2("Diagnostics")
        switch(
            x$index$status,
            absent = cli::cli_alert_info("Access: sequential (no index found)"),
            detected = cli::cli_alert_info("Access: {x$index$type} index detected (not opened)"),
            readable = cli::cli_alert_success("Access: {x$index$type} index readable"),
            unreadable = cli::cli_alert_danger("Access: {x$index$type} index unreadable")
        )

        issues <- x$checks[x$checks$status != "pass", , drop = FALSE]
        if (any(issues$status == "error")) {
            cli::cli_alert_danger("Header checks: failed")
        } else if (nrow(issues)) {
            cli::cli_alert_warning("Header checks: passed with warnings")
        } else {
            cli::cli_alert_success("Header checks: passed")
        }

        if (nrow(issues)) {
            details <- cli::cli_div(class = "inspection-indent")
            for (i in seq_len(nrow(issues))) {
                alert <- switch(
                    issues$status[[i]],
                    warning = cli::cli_alert_warning,
                    error = cli::cli_alert_danger
                )
                alert("{.vcf-id {issues$check[[i]]}}: {issues$message[[i]]}", wrap = TRUE)
            }
            cli::cli_end(details)
        }
    })
}

#' Preview the first few entries in an inspection
#' @param values Character vector of entries.
#' @return One short printable string that includes the remaining count.
#' @noRd
inspection_preview <- function(values) {
    if (!length(values)) { return("none declared") }
    shown <- paste(utils::head(values, 4L), collapse = ", ")
    remaining <- length(values) - 4L
    if (remaining > 0L) paste0(shown, ", ... (+", remaining, ")") else shown
}

#' Create the problem reporter for a single import
#'
#' @param on_problem Reader mode, `"error"` or `"report"`.
#' @return A list of functions. `add` signals a structured error in strict
#'   mode or stores a row in private state for this import. `diagnostics`
#'   returns the diagnostic table.
#' @noRd
new_problem_reporter <- function(on_problem) {
    state <- new.env(parent = emptyenv())
    state$rows <- list(problem_row())
    add <- function(
        chrom,
        pos,
        sample,
        type,
        reason,
        allele_index = NA_integer_,
        allele_offset = NA_integer_,
        value_index = NA_integer_
    ) {
        if (on_problem == "error") {
            locus <- paste0(chrom, ":", if (is.na(pos)) "?" else pos)
            if (!is.na(sample)) { locus <- paste0(locus, " (", sample, ")") }
            cli::cli_abort(
                "{locus}: {reason}",
                class = "methylVCF_import_error",
                chrom = chrom, 
                pos = pos, 
                sample = sample, 
                type = type,
                allele_index = allele_index, 
                allele_offset = allele_offset,
                value_index = value_index,
                reason = reason, 
                call = NULL
            )
        }

        state$rows[[length(state$rows) + 1L]] <- problem_row(
            chrom, 
            pos, 
            sample, 
            type, 
            reason,
            allele_index, 
            allele_offset, 
            value_index
        )

        invisible(NULL)
    }
    list(add = add, diagnostics = function() do.call(rbind, state$rows))
}

#' Construct the structured import diagnostics table
#'
#' @param chrom,pos Chromosome and 1-based position of the source record.
#' @param sample Sample name, or `NA` for a record-level issue.
#' @param type Issue type, either `"invalid"` or `"unsupported"`.
#' @param reason Plain description of the issue.
#' @param allele_index,allele_offset Zero-based allele index and allele offset
#'   when the issue belongs to a particular measurement slot.
#' @param value_index One-based index of the FORMAT measurement, if known.
#' @return A diagnostic data frame. The defaults create its zero-row schema.
#' @noRd
problem_row <- function(
    chrom = character(),
    pos = integer(),
    sample = character(),
    type = character(),
    reason = character(),
    allele_index = integer(),
    allele_offset = integer(),
    value_index = integer()
) {
    data.frame(
        chrom = chrom,
        pos = pos,
        sample = sample,
        type = type,
        allele_index = allele_index,
        allele_offset = allele_offset,
        value_index = value_index,
        reason = reason
    )
}

#' Signal a problem when decoding a VCF record or sample
#'
#' @param reason Plain diagnostic text stored on the condition.
#' @param unsupported Whether the input is well-formed but requires an
#'   unsupported interpretation. Such errors also carry
#'   `methylVCF_unsupported`.
#' @return Signals a `methylVCF_record_error` with the original `reason`.
#' @noRd
abort_methyl_record <- function(reason, unsupported = FALSE) {
    classes <- c(
        if (unsupported) "methylVCF_unsupported", "methylVCF_record_error"
    )
    cli::cli_abort("{reason}", class = classes, reason = reason, call = NULL)
}

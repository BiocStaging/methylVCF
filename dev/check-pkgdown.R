if (!requireNamespace("pkgdown", quietly = TRUE)) {
    cli::cli_abort(
        "Install {.pkg pkgdown} to check the reference index.", call = NULL
    )
}

if (!rmarkdown::pandoc_available()) {
    cli::cli_abort(
        "Install Pandoc or set {.envvar RSTUDIO_PANDOC}.", call = NULL
    )
}

site <- list(destination = tempfile("methylVCF-pkgdown-"))
pkgdown::build_reference_index(override = site)
pkgdown::build_articles_index(override = site)

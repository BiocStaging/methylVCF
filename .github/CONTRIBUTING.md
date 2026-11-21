# Contributing to methylVCF

## Report a bug or suggest a feature

For a bug report, include the methylVCF and caller versions, the genomic
region queried, and a minimal example with the expected and observed
results. Use a small public or synthetic VCF or gVCF that reproduces the
problem.

For a feature request, describe the input, the desired result, and how it
would fit into your workflow.

## Changes and local checks

Edit documentation and `@testexamples` in the roxygen blocks. Regenerate
help and example tests from the package root with `roxygen2::roxygenise()`
and roxytest installed. Include the generated files with their source changes.

Keep regressions in the `tests/testthat` file for the behavior they cover.
Use `devtools::test()` to run the suite locally. Use the small fixtures in
`inst/extdata` for examples and tests, and update the corresponding
[generation recipe][fixtures] when a fixture changes. The separate caller
validation script in `dev/` uses externally obtained VCFs.

Edit `README.Rmd` and regenerate `README.md` from the package root with
rmarkdown, pkgload, Pandoc and the package dependencies installed.
Include both files with README changes.

```sh
Rscript --vanilla -e 'rmarkdown::render("README.Rmd")'
```

Run `R CMD build .`, then `R CMD check` and `BiocCheck::BiocCheck()` on the
resulting tarball. When adding a public topic or vignette, update
`_pkgdown.yml` and run the reference check with pkgdown and Pandoc installed.

```sh
Rscript --vanilla dev/check-pkgdown.R
```

## Open a pull request

Open the pull request against `devel`, link any related issue, and describe
what changed and which checks you ran. GitHub Actions runs Bioconductor
devel checks on Linux, macOS, and Windows, including TBI and CSI queries
through Rsamtools.

[fixtures]: https://github.com/xec-cm/methylVCF/blob/devel/inst/scripts/README.md

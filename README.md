
<!-- README.md is generated from README.Rmd. Edit README.Rmd. -->

# methylVCF <img src="man/figures/logo.png" align="right" width="140" alt="" />

<!-- badges: start -->

[![R CMD check on
devel](https://github.com/xec-cm/methylVCF/actions/workflows/check-bioc.yml/badge.svg?branch=devel)](https://github.com/xec-cm/methylVCF/actions/workflows/check-bioc.yml?query=branch%3Adevel)
[![MIT
license](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/license/mit)
<!-- badges: end -->

**DNA methylation from VCF and gVCF, by genomic region.**

methylVCF reads 5mC measurements reported by Illumina DRAGEN 5-Base and
Rastair. Work with individual measurements as `GRanges`, or assemble
reference methylation across samples in a `RangedSummarizedExperiment`.

## Installation

Requires R 4.6 or later.

### Development version

Install from [GitHub](https://github.com/xec-cm/methylVCF) with pak.

``` r
pak::pak("xec-cm/methylVCF")
```

### Bioconductor

This installation will be available once methylVCF is published on
Bioconductor.

``` r
BiocManager::install("methylVCF")
```

## Quick start

The examples use plyranges for genomic ranges and
tidySummarizedExperiment for site-by-sample tables. Install the
companions to run the tutorials.

``` r
BiocManager::install(c("plyranges", "tidySummarizedExperiment"))
```

Use `inspectMethylVcf()` to check sample names, methylation fields and
index readability before a regional query.

Build a one-base interval with plyranges and read the reported
measurements from a bundled DRAGEN gVCF.

``` r
library(methylVCF)

vcf <- system.file("extdata", "dragen-multisample.g.vcf", package = "methylVCF")
site <-
    tibble::tibble(seqnames = "chr21", start = 13204474L, width = 1L) |>
    plyranges::as_granges()

readMethylVcf(vcf, site) |>
    tibble::as_tibble() |>
    dplyr::select(seqnames, start, sample, fraction, informative_depth)
#> # A tibble: 2 × 5
#>   seqnames    start sample   fraction informative_depth
#>   <fct>       <int> <chr>       <dbl>             <int>
#> 1 chr21    13204474 sample_a    0.238                21
#> 2 chr21    13204474 sample_b    0                     0
```

Use `plyranges::join_overlap_inner()` to annotate calls with genomic
features. For aligned reference measurements, pass the original ranges
to `asMethylExperiment()`. The experiment keeps reported values and
coverage-backed fractions in separate assays.

## Articles

- **Explore methylation in a region** with `vignette("methylVCF")`.
- **Inspect and read your own DRAGEN or Rastair files** with
  `inspectMethylVcf()` and `vignette("input-vcf")`.
- **Reproduce the HG002 case study** with `vignette("dragen-hg002")`.

## Community

- [Report a
  bug](https://github.com/xec-cm/methylVCF/issues/new?template=bug-report.md).
  Include the caller version and a small reproducible genomic query.
- [Suggest a
  feature](https://github.com/xec-cm/methylVCF/issues/new?template=feature_request.md).
  Describe the input and what you want to obtain from it.
- [Contribute to methylVCF](.github/CONTRIBUTING.md). See how to work on
  the package and run the checks locally.

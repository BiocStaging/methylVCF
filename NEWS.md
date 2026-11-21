# methylVCF 0.99.0

- Inspect VCF headers with `inspectMethylVcf()` to check sample IDs,
  declared fields and index readability before making a regional query.
- Query reported 5mC by genomic region and sample with `readMethylVcf()`.
  Supports Illumina DRAGEN 5-Base VCF and gVCF, including reference blocks
  and variants, and Rastair VCF with identifiable reference and de novo
  CpG measurement slots.
- Return one-base measurement anchors as `GRanges` with reported fractions,
  read depths, coverage state, source record coordinates and genotype
  evidence. Preserve DRAGEN masks and assign alleles and phased slots
  where supported by the record.
- Query local BGZF inputs through Rsamtools with TBI or CSI indexes.
  Retrieve disjoint regions directly and scan unindexed inputs sequentially.
  Preserve the VCF sequence dictionary and declared assemblies in `Seqinfo`.
- Expose import diagnostics for malformed or unsupported measurements.
  Report mode retains readable rows and records issues such as
  unprojectable allele offsets, ambiguous Rastair slot layouts, and slots
  with more modified reads than informative reads.
- Use `asMethylExperiment()` to align reference sites across samples in a
  `RangedSummarizedExperiment`. Keep reported fractions, coverage-backed
  fractions, depth and coverage status separately, with sample and import
  provenance.
- Include offline examples for both callers and a reproducible HG002 DRAGEN
  case study.

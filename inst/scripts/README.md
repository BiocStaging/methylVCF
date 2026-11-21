# Reproduce the bundled VCF examples

The files under [extdata](../extdata/README.md) are small, local inputs for
methylVCF tests and tutorials. These recipes run in a scratch directory.
The packaged examples are available without running a caller or downloading
sequencing data.

## Curated VCF records

Run these commands from the root of a checkout. After installing methylVCF,
you can locate the same script in R with
`system.file("scripts", "generate-vcf-examples.py", package = "methylVCF")`.

```sh
package_root=$(pwd)
scratch=$(mktemp -d)
cd "$scratch"
python3 "$package_root/inst/scripts/generate-vcf-examples.py"
```

- `dragen-example.g.vcf` reconstructs the block in the
  [Illumina methylation annotation example][dragen]. The REF mask marks an
  anchor at POS + 3, where the fraction is 0.2381 and the depth is 21.
  The artificial variants cover both phased genotype orders, an absent
  fraction, and pooled and phased homozygous measurements.
- `dragen-multisample.g.vcf` reuses that reference block and adds two
  artificial samples, an observed zero, a zero-depth report, and an absent
  sample measurement. The two samples have reversed genotype order at the
  same phased variant.
- `rastair-example.vcf` was written against the Rastair 2.3.0-rc1
  [VCF header snapshot][header] and [per-context regression test][contexts].
  It has forward and reverse reference CpG measurements, a site with two
  distinct context slots, and a de novo context with ambiguous local
  allele assignment.

From the scratch directory, use `cmp "$package_root/inst/extdata/FILE" FILE`
with each generated file name to check byte-for-byte equivalence.

## VCF emitted by Rastair 2.2.0

Use Python 3, samtools and Rastair 2.2.0. Continue in the same scratch
directory. The generator uses a fixed seed to construct an artificial
reference and single-end alignments. Rastair generates the calls from
those inputs.

```sh
python3 "$package_root/inst/scripts/generate-rastair-fixture.py"
samtools faidx synthetic.fa
samtools view -b -o synthetic.bam synthetic.sam
samtools index synthetic.bam
rastair call --threads 2 --fasta-file synthetic.fa --unpaired --no-ml --all \
    --vcf-info-fields SC5,CPG,CPGnovo \
    --vcf-format-fields GT,M5mC,DPM5mC,ADM5mC \
    --vcf synthetic.vcf synthetic.bam
cmp "$package_root/inst/extdata/rastair-2.2.0-synthetic.vcf" synthetic.vcf
```

This output includes a record at position 302 with both CpG flags and
one measurement slot, which exercises the unsupported-record diagnostic.
The packaged file has SHA-256
`2f78e9427dfbf8e3d6262b05b9a79930a73f38cb3a121b5a6a64dda1225d40f7`.
Rastair software attribution is Copyright © 2025, University of Oxford.
[Rastair license][license]

The [caller validation report][validation] covers the separate HG002 and
biological Rastair comparisons, including the input sources and expected
results.

[dragen]: https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/core-functionality/methylation-annotation/
[header]: https://github.com/bsbludwig/rastair/blob/23d096617e3f9f46bd9405b16518b2111868c02b/tests/snapshots/call_cli__vcf_with_all_fields.snap
[contexts]: https://github.com/bsbludwig/rastair/blob/23d096617e3f9f46bd9405b16518b2111868c02b/tests/call_records.rs#L73-L110
[license]: https://github.com/bsbludwig/rastair/blob/23d096617e3f9f46bd9405b16518b2111868c02b/LICENSE.txt
[validation]: https://github.com/xec-cm/methylVCF/blob/devel/dev/caller-validation.md

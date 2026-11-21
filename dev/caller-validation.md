# Validation against caller output

This check compares VCF fields with `readMethylVcf()` and, for reference
CpGs, with all five assays in `asMethylExperiment()`. The package tests and
vignette use small offline fixtures.

## Rastair 2.2.0 on mouse TAPS reads

The [Rastair upstream fixtures][inputs] contain an indexed BAM and BGZF
reference. The BAM identifies sample `mES_TAPS` and records alignment to
mm10 with spike-ins. The `chr19` reference is 61,431,566 bases long in both
the BAM sequence dictionary and the FASTA index. The VCF sample emitted by
Rastair is named `sample`.

Use a binary for release 2.2.0. The check used the [Bioconda Apple Silicon
build][bioconda] with SHA-256
`dec9f37cd7c19e478f8fabec7d082250c56c15187c9de25e23b5f60e42a342c3`.
Keep the BAM, FASTA and their companion index files in the working directory.
Run Rastair with the following arguments.

```sh
rastair call --threads 2 --fasta-file test.fasta.gz \
    --region chr19:6103156-6143156 --no-ml --all --cpgs-only \
    --vcf-info-fields SC5,CPG,CPGnovo \
    --vcf-format-fields GT,M5mC,DPM5mC,ADM5mC \
    --vcf biological.vcf test.bam

bcftools reheader -f test.fasta.gz.fai \
    -o biological-reheader.vcf biological.vcf
```

The regional output from 2.2.0 declares `chr19` length 40,000 even though it
writes chromosome-relative positions. The reheader step replaces the
declaration using the aligned FASTA index. All 1,901 VCF data lines were
byte-identical before and after this step. The reader preserves the supplied
sequence dictionary. Without the correction, `GRanges` reports that these
positions exceed the declared length.

Validate the output from the root of an installed checkout. The script
takes the reheadered Rastair VCF and the indexed public HG002 gVCF or the
local excerpt described in the vignette.

```sh
Rscript --vanilla dev/validate-caller-output.R \
    /path/to/biological-reheader.vcf /path/to/HG002.gvcf.gz
```

The raw FORMAT values matched all **1,901** measurements, comprising 949
forward-strand and 952 reverse-strand slots. All were reference CpG contexts
with allele index 0 and missing haplotype. The result retained 457 filtered
records. The experiment kept **1,444** PASS sites and recovered the same
site order, fractions, depths, coverage states and VCF sample ID as the
source. At chr19:6106281, the caller reports a zero fraction with zero
informative depth. The experiment keeps the raw 0 and sets the analytical
fraction to missing.

| POS | REF | ALT | FILTER | M5mC | DPM5mC | ADM5mC |
| --- | --- | --- | --- | ---: | ---: | ---: |
| 6105744 | C | T | PASS | 0.8 | 10 | 8 |
| 6105745 | G | . | PASS | 0.5 | 6 | 3 |
| 6106281 | C | . | PASS | 0 | 0 | 0 |

The default validation uses the CpG-only output. A separate `--all` run
without `--cpgs-only` emitted 1,080 records with methylation values and
neither CpG flag in this interval. Those records lack the context required
by the importer. The package's interpretation covers flagged CpG records
with an identifiable slot layout. [Rastair VCF format][format]

| File | SHA-256 |
| --- | --- |
| `test.bam` | `90d52f3cf8a119ddd0b3cabfbc4e98e231e98847fb12cc12bdf130e5fa0a60be` |
| `test.bam.bai` | `2553ad5a199ca6a49f0fecc4be39ad93b122c07bd88f6cc58d8571f913c7166e` |
| `test.fasta.gz` | `626d36fffc555deee325c10fd417d9b3b8f5cbed6b8f9a3f1cad6af7d0788259` |
| `test.fasta.gz.fai` | `fbae63c6f2364e3d0e7ba8c6f6d7e2e76c1d5315c3ea14c5da8d2dedfdeec380` |
| `test.fasta.gz.gzi` | `bd2112fb219df2657857fb4eca63c8b849d13e07f6e994f68ef4e0b626620a10` |
| `biological.vcf` | `6487aa0f45543f41ed7603fcadf4740520e5702538427d4d6d1ea09c7db0c042` |
| `biological-reheader.vcf` | `6fd60dc0c80eedace6ec2690aed3b10eb720c222adf08cb57a60cec0a098e87b` |

## Rastair 2.2.0 on artificial reads

The packaged VCF was generated on artificial reads. The
[fixture recipe](../inst/scripts/README.md) records the generator, command
and checksum.

At positions 101 and 102, the two strands of a reference CpG report
0.666667 with depth 12 and eight modified reads. At position 301, a
neighboring substitution forms a de novo CpG anchored on the local REF C.
At position 302, a biallelic C-to-G record has both `CPG` and `CPGnovo`,
and each methylation field contains one scalar. The record can have two
different measurements. Reporting it as unsupported avoids giving its
value to either context. Position 303 reports zero fraction at zero depth.
Separate hand-written fixtures cover explicit missing fractions and both
slots at one position.

The fixture uses synthetic sequence and reads created for methylVCF.
The upstream biological files are fetched only for local comparison.
The upstream repository makes them publicly accessible without a separate
redistribution license identified in its data directory. Rastair software
attribution is Copyright © 2025, University of Oxford. [Rastair license][license]

## DRAGEN 5-Base on HG002

The public [Illumina BaseSpace demo][hg002] has one germline sample
`NA24385` on GRCh38 from DRAGEN 4.4.6. The indexed local excerpt used here
has SHA-256 `1ab9ccee85893dfeaedcfef4c8bf99eb47d5e9f5639f080cb2867bd4ed43ab86`.
The validation script scans the raw VCF lines and checks the masks and
FORMAT vectors independently of the package's indexed query.

| chr1 anchor | Record POS | GT | REF/ALT index | Haplotype | M5mC | DPM5mC |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 10493 | 10486 | 0/0 | 0 | NA | 0.4839 | 31 |
| 10816 | 10816 | 0/0 | 0 | NA | 0 | 0 |
| 28692 | 28692 | 0&#124;1 | 0 | 1 | 0 | 14 |
| 94011 | 94011 | 0/1 | 1 | NA | 0.9286 | 14 |
| 171535 | 171535 | 1&#124;1 | 1 | NA | 0.7692 | 13 |

At chr1:3042652, the source genotype is `0/3`. The caller reports one REF
slot and four ALT slots. The reader retains the two measurements at the
shared anchor and reports unprojectable ALT 3 offsets 5, 19 and 31 in
FORMAT slots 3, 4 and 5. [DRAGEN methylation masks][dragen]

[inputs]: https://github.com/bsbludwig/rastair/tree/23d096617e3f9f46bd9405b16518b2111868c02b/tests/data
[bioconda]: https://anaconda.org/bioconda/rastair/files
[format]: https://www.rastair.com/formats/vcf.html
[license]: https://github.com/bsbludwig/rastair/blob/23d096617e3f9f46bd9405b16518b2111868c02b/LICENSE.txt
[hg002]: https://basespace.illumina.com/analyses/864531910/files/ds.48b7596730dd47ef97699b59ccd3641d
[dragen]: https://illumina.github.io/IlluminaConnectedAnnotationsDocumentation/core-functionality/methylation-annotation/

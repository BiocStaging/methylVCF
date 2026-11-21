"""Write the three curated VCF examples into the current directory."""

from pathlib import Path


def write_vcf(path, metadata, samples, records):
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    lines = metadata.strip().splitlines()
    lines.append("\t".join(columns + samples))
    lines.extend("\t".join(record.split()) for record in records)
    Path(path).write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


dragen_metadata = """
##fileformat=VCFv4.5
##source=DRAGEN
##contig=<ID=chr21,length=46709983>
##ALT=<ID=NON_REF,Description="Represents any unobserved alternate allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="Last reference position of a block">
##INFO=<ID=M5mC,Number=R,Type=String,Description="Per-allele mask of reportable modification positions">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
##FORMAT=<ID=M5mC,Number=M,Type=Float,Description="5mC fraction">
##FORMAT=<ID=DPM5mC,Number=M,Type=Integer,Description="Informative 5mC depth">
"""

# A 42-base reference block with the reported CpG at POS + 3.
mask = "...Z" + "." * 38
block = (
    f"chr21 13204471 . A <NON_REF> . PASS END=13204512;M5mC={mask},. "
    "GT:M5mC:DPM5mC 0/0:0.2381:21"
)
phased_snv = (
    "chr21 13204520 . C G,<NON_REF> 60 PASS M5mC=x,x,. "
    "GT:M5mC:DPM5mC:PS 1|0:0.8,0:9,14:13204520"
)

write_vcf(
    "dragen-example.g.vcf",
    dragen_metadata,
    ["dragen_demo"],
    [
        block,
        phased_snv,
        "chr21 13204534 . C G,<NON_REF> 60 PASS M5mC=x,x,. "
        "GT:M5mC:DPM5mC:PS 0|1:.,0.7:0,10:13204520",
        "chr21 13204540 . T C,<NON_REF> 60 PASS M5mC=.,Z,. "
        "GT:M5mC:DPM5mC:PS 1|1:0.4:15:13204520",
        "chr21 13204550 . T C,<NON_REF> 60 PASS M5mC=.,Z,. "
        "GT:M5mC:DPM5mC:PS 1|1:0.2,0.6:5,10:13204520",
    ],
)
write_vcf(
    "dragen-multisample.g.vcf",
    dragen_metadata,
    ["sample_a", "sample_b"],
    [
        block + " 0/0:0:0",
        "chr21 13204515 . C <NON_REF> . PASS M5mC=Z,. "
        "GT:M5mC:DPM5mC 0/0:0:9 ./.:.:.",
        phased_snv + " 0|1:0.3,0.7:10,10:13204520",
    ],
)

rastair_metadata = """
##fileformat=VCFv4.5
##rastairVersion=2.3.0-rc1
##contig=<ID=chr19,length=61431566>
##INFO=<ID=SC5,Number=1,Type=String,Description="Five-base reference context centered on POS">
##INFO=<ID=CPG,Number=0,Type=Flag,Description="Reference CpG context">
##INFO=<ID=CPGnovo,Number=0,Type=Flag,Description="Candidate de novo CpG context">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=M5mC,Number=.,Type=Float,Description="5mC fraction in CpG context order">
##FORMAT=<ID=DPM5mC,Number=.,Type=Integer,Description="Informative 5mC depth in CpG context order">
##FORMAT=<ID=ADM5mC,Number=.,Type=Integer,Description="5mC-supporting depth in CpG context order">
"""

write_vcf(
    "rastair-example.vcf",
    rastair_metadata,
    ["rastair_demo"],
    [
        "chr19 6105744 . C . . PASS CPG;SC5=AACGT "
        "GT:M5mC:DPM5mC:ADM5mC 0/0:0.8:10:8",
        "chr19 6105745 . G . . PASS CPG;SC5=ACGTA "
        "GT:M5mC:DPM5mC:ADM5mC 0/0:0.5:14:7",
        "chr19 6105750 . C G 60 PASS CPG;CPGnovo;SC5=ACCGC "
        "GT:M5mC:DPM5mC:ADM5mC 1/1:.,1:0,5:0,5",
        "chr19 6105760 . C G 60 PASS CPGnovo "
        "GT:M5mC:DPM5mC:ADM5mC 0/1:0.6:5:3",
    ],
)

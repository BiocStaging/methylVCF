from pathlib import Path
import random

random.seed(817)
sequence = list(random.choices("AT", k=700))
for site, motif in ((101, "CG"), (301, "CCG"), (501, "CG")):
    sequence[site - 1:site + len(motif) - 1] = motif
reference = "".join(sequence)
Path("synthetic.fa").write_text(f">synthetic\n{reference}\n")

# Aligned single-end reads on both strands. The middle C in CCG is G in
# the second window, where forward/reverse reads also show the conversion.
with Path("synthetic.sam").open("w") as sam:
    print("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:synthetic\tLN:700",
          "@RG\tID:artificial\tSM:artificial", sep="\n", file=sam)
    for target in (101, 302):
        for offset in range(12):
            start = target - 48 + offset
            seq = reference[start - 1:start + 99]
            base = target - start
            for flag in (0, 16):
                read = list(seq)
                if target == 101 and offset >= 4:
                    read[base if flag == 0 else base + 1] = "T" if flag == 0 else "A"
                if target == 302:
                    read[base] = "G"
                    if flag == 0 and offset >= 3:
                        read[base - 1] = "T"
                name = f"r{target}_{offset}_{flag}"
                print(name, flag, "synthetic", start, 60, "100M", "*", 0, 0,
                      "".join(read), "I" * 100, "RG:Z:artificial", sep="\t", file=sam)

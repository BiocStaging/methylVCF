# Indexed VCF fixtures

These BGZF files contain the same synthetic records as the plain VCFs in
`inst/extdata/`. Each file has a CSI companion. The tests build TBI indexes
with Rsamtools and use these fixtures to exercise CSI.

To regenerate the CSI fixtures from the package root, use HTSlib's
`tabix -C -p vcf` after compressing the corresponding input with
`Rsamtools::bgzip()`.

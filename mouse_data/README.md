# Mouse reference panel

## Existing RDS objects

No existing RDS file is changed by the PLINK build.

| File | Class and dimensions | IDs and contents |
|---|---|---|
| `genotypes.rds` | `data.frame`, 1,267 × 1,813 | Rows and columns are uniquely named `1`–`1267` and `1`–`1813`. Numeric hard calls are exactly 0/1/2 with 19,624 missing values. Per-individual missing counts range from 0 to 1,813 and per-marker counts from 8 to 55. |
| `genotypes_imputed.rds` | numeric `matrix`, 1,267 × 1,813 | Same IDs and order. All 2,297,071 entries are exact integer-valued 0/1/2 calls; there are no continuous dosages and no missing values. All observed calls in `genotypes.rds` are unchanged. The formerly missing calls comprise 5,585 zeros, 11,152 ones, and 2,887 twos. |
| `mouse.rds` | `data.frame`, 1,177 × 6 | Row IDs `91`–`1267`; columns `sire`, `dam`, `sex`, `reps`, `Gl`, `BW`; no missing values. `sire`, `dam`, `sex`, and `reps` are factors; `Gl` and `BW` are numeric phenotypes. |
| `mouseqtl.rds` | `data.frame`, 1,177 × 8 | The six `mouse.rds` columns plus factor genotypes `M227` and `M1139`, coded `AA`/`AB`/`BB`; 3 and 13 missing calls respectively. |
| `pedigree.rds` | `data.frame`, 1,267 × 6 | Row order and integer `id` exactly match genotype rows. Columns are `id`, `sire`, `dam`, `family`, `sex`, `generation`; no missing values. Parent 0 denotes unknown/founder. Sex is `Male`/`Female`; generations are `M6`, `IC`, `F1`, and `F2`; `family` has 68 observed sire/dam labels. |

The RDS objects contain no marker chromosome, base-pair, genetic-map, or
nucleotide allele metadata. `mouse.rds` and `mouseqtl.rds` provide phenotypes
for only 1,177 animals, and the repository does not identify either `Gl` or
`BW` as the canonical PLINK phenotype. The FAM phenotype is therefore `-9`.

## Canonical PLINK files

`mouse.bed`, `mouse.bim`, and `mouse.fam` use
`genotypes_imputed.rds`, because it is a fully observed hard-call matrix—not a
continuous dosage matrix—and preserves every observed call from
`genotypes.rds`. No rounding occurs.

Individual and marker order are unchanged. FAM family, paternal, maternal, and
sex fields come from `pedigree.rds`; male is PLINK code 1 and female code 2.
Because marker biology is absent, BIM chromosome, genetic map, and base-pair
position are explicit unknown values (`0`). Alleles `1` and `2` are
non-biological coding labels. Marker IDs are the original unique numeric column
names. Do not interpret any placeholder as a coordinate or allele.

The builder uses `genio` 1.1.2 and fails on non-integer dosages. It validates
the temporary fileset through both `genio` and `qgg::gprep()`/`qgg::getG()`
before installing it. Exact checksums, package versions, and results are in
`mouse_validation.txt` and `mouse_sessionInfo.txt`.

```text
Rscript scripts/build_mouse_plink.R
Rscript scripts/build_mouse_plink.R --force  # intentional rebuild only
```

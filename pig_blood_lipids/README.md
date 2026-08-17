# Pig blood lipids

This dataset contains real GWAS genotypes and six blood-lipid traits from
Laiwu, Erhualian, and Duroc x (Landrace x Yorkshire; DLY) pigs. Populations are
deliberately retained as separate PLINK datasets; no cross-population marker or
allele harmonization is assumed.

## Contents

| Population | PLINK prefix | Individuals | Markers |
|---|---|---:|---:|
| Laiwu | `laiwu/laiwu` | 316 | 61,565 |
| Erhualian | `erhualian/erhualian` | 334 | 61,565 |
| Duroc x (Landrace x Yorkshire) | `dly/dly` | 610 | 61,565 |

`phenotypes.csv` contains FID, IID, population, TCHOL, TG, HDL-C, LDL-C,
HDL-C/LDL-C, and AI. Archived values are retained as text and missing values as
`NA`; they are not imputed, filtered, standardized, adjusted, or transformed.
Per-population nonmissing counts and numerical summaries are in
`validation.txt`.

The paper reports 1,256 pigs analyzed, while the archived GWAS PED files contain
1,260 genotype rows (316 Laiwu, 334 Erhualian, and 610 DLY). The canonical
phenotype table follows genotype IDs exactly: all Laiwu IDs match, one
Erhualian genotype (`E1481B`) has no phenotype row and therefore has six `NA`
traits, and 89 DLY phenotype-only rows are not installed because they have no
corresponding genotype. No genotype row was removed.

PLINK2 converts each archived PED/MAP pair without marker QC. FID, IID, parents,
sex, phenotype field, missing calls, marker IDs, positions, and marker order are
preserved. The archive uses numeric PLINK aliases 23/24 for X/Y; PLINK2 writes
the canonical `X`/`Y` BIM labels. An absent allele is written as `.` rather than
the source PED missing code `0`. Full streaming PED round trips verify each
unordered allele pair. The archive does not identify a genome assembly, so no
build is assigned here.

## Provenance and licence

Source: Yang, Hui; Huang, Xiaochang; Zeng, Zhijun; Zhang, Wanchang; Liu,
Chenlong; Fang, Shaoming; Huang, Lusheng; Chen, Congying (2016). *Data from:
Genome-wide association analysis for blood lipid traits measured in three pig
populations revealed a substantial level of genetic heterogeneity* [Dataset].
Dryad. <https://doi.org/10.5061/dryad.4gh70>

Associated paper: Yang et al. (2015), *PLOS ONE*,
<https://doi.org/10.1371/journal.pone.0131667>.

Dryad version endpoint:
<https://datadryad.org/api/v2/versions/25357/download>. The inner archive is
`Primary_data.zip` (22,332,040 bytes; MD5
`51bcd5b1011b55379a0425bf5c0e2a2c`, SHA-256
`71b2014f71c9fc2cbd2f1923fa1da9c46d2d9523041292ef7a4e6966f6e50a05`).
Dryad releases the dataset under CC0 1.0.

Build and validate:

```text
Rscript scripts/build_livestock_data.R --species=pig
Rscript scripts/validate_livestock_data.R
```

Existing output is never replaced unless `--force` is supplied to the build.

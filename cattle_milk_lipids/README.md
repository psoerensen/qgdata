# Cattle milk polar lipids

This dataset contains real high-density imputed genotypes and milk
polar-lipid measurements for 336 dairy cows. The archived PLINK files are
copied byte-for-byte and only renamed to the canonical `cattle_milk_lipids`
basename.

## Contents

- `cattle_milk_lipids.{bed,bim,fam}`: 336 individuals and 593,870 markers;
  SNP-major PLINK1 binary data, with no missing genotype calls.
- `phenotypes.csv`: the unchanged source CSV (FID, IID, and 56 traits).
- `fixed_effects.csv`: the source FID, IID, and combined `Year_Batch` effect,
  represented as CSV; its nine archived levels are B1--B9.
- `validation.txt`: dimensions, missingness, trait summaries, checksums,
  software versions, and focused qgg/genio round trips.

Trait prefixes identify the source groups PS, SM, PE, PC, PI, LPC, LaCer, and
Glucer. The source documentation expands PS, SM, PE, PC, PI, LaCer, and Glucer
as phosphatidylserine, sphingomyelin, phosphatidylethanolamine,
phosphatidylcholine, phosphatidylinositol, lactosylceramide, and
glucosylceramide. It does not define measurement units or expand LPC, so this
repository does not infer them. The source states that outliers have not been
removed. Missing phenotype cells are retained as `NA`; per-trait counts and
ranges are recorded in `validation.txt`.

The source describes the genotypes as HD imputed data based on Run7 of the
1000 Bull Genomes Project. The archive does not state a genome assembly, so no
genome build is assigned here. Biological BIM chromosomes, positions, and
alleles are retained exactly. These data are useful for quantitative-genetics
examples involving real LD and multiple correlated lipid traits, but the small
sample size and undocumented measurement units/build should be considered.

## Provenance and licence

Source: Ghoreishifar, Mohammad; Macleod, Iona; Chamberlain, Amanda; Liu,
Zhiqian; Lopdell, Thomas; Littlejohn, Mathew; Xiang, Ruidong; Pryce, Jennie;
Goddard, Michael (2025). *Data from: An integrative approach to prioritize
candidate causal genes for complex traits in cattle* [Dataset]. Dryad.
<https://doi.org/10.5061/dryad.bcc2fqzph>

Dryad version endpoint:
<https://datadryad.org/api/v2/versions/368250/download>. The inner archive is
`Polar_Lipids_Phenotypes_HDGenotypes.zip` (16,745,878 bytes; SHA-256
`d92af1f1ddd63b401b18e88d6b65ed20cdc5c366cdff24fde16224c1d34f9590`).
Dryad releases the dataset under CC0 1.0.

Build and validate:

```text
Rscript scripts/build_livestock_data.R --species=cattle
Rscript scripts/validate_livestock_data.R
```

Existing output is never replaced unless `--force` is supplied to the build.

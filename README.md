# qgdata

`qgdata` is a small collection of fixed reference genotype panels for
`sblrbench`, `sblr`, `qgg`, and related quantitative-genetics work. Phenotypes
and genetic architectures are simulated elsewhere with `gsim()`.

## Datasets

| Dataset | Individuals x markers | Status |
|---|---:|---|
| `mouse_data/mouse` | 1,267 x 1,813 | Available and validated. |
| `human_1000g_eur/human_1000g_eur` | 5,000 x 50,000 | Available and validated; EUR mosaic rate 20. |
| `human_1000g_afr/human_1000g_afr` | 5,000 x 50,000 | Available and validated; AFR mosaic rate 25. |
| `cattle_milk_lipids/cattle_milk_lipids` | 336 x 593,870 | Real HD imputed cattle genotypes with 56 milk polar-lipid traits. |
| `pig_blood_lipids/laiwu/laiwu` | 316 x 61,565 | Real Laiwu genotypes with six blood-lipid traits. |
| `pig_blood_lipids/erhualian/erhualian` | 334 x 61,565 | Real Erhualian genotypes with six blood-lipid traits. |
| `pig_blood_lipids/dly/dly` | 610 x 61,565 | Real Duroc x (Landrace x Yorkshire) genotypes with six blood-lipid traits. |
| `simulated_human_data/human` | 5,000 x 50,000 | Protected legacy independent-marker fixture. |
| `1000G/sample_chr1` | 489 x 1,000 | Protected legacy chromosome 21 fixture. |
| `1000G/sample_chr2` | 489 x 1,000 | Protected legacy chromosome 22 fixture. |

The two human ancestry panels are synthetic reference-conditioned genotypes,
not original 1000 Genomes individuals and not a demographic or privacy model.
Each synthetic haplotype copies phased alleles from same-ancestry donors and
switches donors between adjacent markers with probability
`1 - exp(-copy_rate * delta_morgan)`. There is no mutation, ancestry mixing,
frequency adjustment, imputation, artificial LD blocking, or post-generation
individual filtering.

## Human-panel reproduction and validation

The locked production command is:

```text
Rscript scripts/build_human_mosaic.R --mode=full --ancestry=both
```

It uses the final 20130502 Phase 3 GRCh37 chromosome 22 phased reference, the
deterministic shared 50,000-marker design, EUR rate 20, AFR rate 25, and seeds
20960818 and 20960819. Blocks are checkpointed under ignored
`cache/human_r_mosaic/production/`; use `--resume` after interruption. `--force`
explicitly discards compatible build checkpoints and is not equivalent to
`--resume`. Canonical files are installed only after both staged panels pass.

Run the maintained focused checks with:

```text
Rscript scripts/test_haplotype_mosaic.R
Rscript scripts/validate_human_panels.R
```

Production validation uses full-panel frequency, duplicate, hard-call, and
PLINK2 KING summaries plus deterministic 10,000-marker scientific subsets. It
does not form a dense 50,000 x 50,000 LD matrix. Five synthetic 200-person
subsets are compared with 20 deterministic matched 200-person reference
subsets. AFR q99.5 kinship is retained as a descriptive upper-tail warning;
there were zero duplicate/MZ and first-degree pairs, so it is not an isolated
hard failure. See `validation/human_validation_report.md` and the compact CSV
tables beside it.

Historical generator comparison is intentionally short: `sim1000G` preserved
MAF but failed local-LD fidelity in the smoke study; HAPNEST was prepared but
not evaluated because its execution environment and custom-reference inputs
were unavailable. Neither superseded execution workflow is maintained.

## Source and citation

The source is the official phased 1000 Genomes Project Phase 3 final release,
GRCh37 chromosome 22:

- VCF: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz>
- sample panel: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel>
- Phase 3: <https://www.internationalgenome.org/data-portal/data-collections/phase3/>
- IGSR data use: <https://www.internationalgenome.org/faq/do-i-need-permission-to-use-igsr-data-in-my-own-scientific-research/>

EUR donors are CEU, GBR, FIN, IBS, and TSI. AFR donors are ESN, GWD, LWK,
MSL, and YRI; ASW and ACB are excluded. Cite the 1000 Genomes Project
Consortium, *Nature* 526, 68-74 (2015),
<https://doi.org/10.1038/nature15393>, and Clarke et al., *Nucleic Acids
Research* 45, D854-D859 (2017), <https://doi.org/10.1093/nar/gkw829>.

Downloaded sources, caches, checkpoints, temporary files, and candidate panels
are ignored. Intended canonical BED/BIM/FAM files are not ignored.

## Livestock datasets

`cattle_milk_lipids` contains the source authors' unchanged PLINK files for 336
cows, 56 milk polar-lipid phenotypes, and the combined year/batch fixed effect
(nine levels). The genotypes are real high-density imputed genotypes based on
Run7 of the 1000 Bull Genomes Project. See
[`cattle_milk_lipids/README.md`](cattle_milk_lipids/README.md).

`pig_blood_lipids` contains the archived GWAS data for Laiwu, Erhualian, and
Duroc x (Landrace x Yorkshire) pigs. The three populations remain separate
because their marker/allele representations must not be assumed harmonized.
The shared phenotype table contains FID, IID, population, and the six archived
blood-lipid traits. See [`pig_blood_lipids/README.md`](pig_blood_lipids/README.md).

Reproduce and validate with:

```text
Rscript scripts/build_livestock_data.R --species=both
Rscript scripts/validate_livestock_data.R
```

The builder downloads the anonymous Dryad version archives into ignored
`cache/livestock/`, follows redirects, verifies the pinned inner-archive
checksums, and requires `--force` before replacing existing outputs. Both
datasets are released by Dryad under CC0; their README files give the full
recommended citations and DOI links.

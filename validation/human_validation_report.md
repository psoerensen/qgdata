# Human ancestry-panel production validation

Status: both panels are available and validated as of 2026-08-17.

## Design and provenance

The panels use phased chromosome 22 haplotypes from the final 20130502 1000
Genomes Phase 3 release on GRCh37. EUR donors are CEU, GBR, FIN, IBS, and TSI;
AFR donors are ESN, GWD, LWK, MSL, and YRI, excluding ASW and ACB. The same
50,000 biallelic SNPs, IDs, alleles, positions, and order are used for both.

Each haplotype starts from a uniformly sampled same-ancestry donor. Between
markers it switches to a different uniformly sampled donor with probability
`1 - exp(-copy_rate * delta_morgan)`. Independently generated haplotypes are
paired into 0/1/2 diploid hard calls. EUR uses rate 20 and seed 20960818; AFR
uses rate 25 and seed 20960819. There is no mutation, ancestry mixing,
frequency adjustment, imputation, artificial blocking, or individual filter.

Original VCF SHA-256 is
`a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5`.
The cached phased PGEN SHA-256 is
`9e07c2f699d21a449917377e8ac6831aacdb9bddfbb8431822a9ab164a9caac1`;
the PVAR SHA-256 is
`0661d9202fbd71c88439e73cefd94a3cab6e77a24c737e014b7f6cdb581a8af5`;
the PSAM SHA-256 is
`909387e1783b6f759667430e003fadcd439c4ea18c148e0f81a73f15a7478141`;
the population-panel SHA-256 is
`b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858`;
the map SHA-256 is
`057f35f9e53bcaca83b04fb9a4c76948f19b0f5cef0d88b548314b9670245d50`;
and the selected-variant file SHA-256 is
`d8cc5e061c3ed2a6b9bbf06eb057833b9f96eabf846f7c1d7f53431b24b217e6`.

## Scaling decision

The locked 1,000 x 10,000 checkpoint was accepted for both ancestries without
rerunning or retuning. EUR rate 20 passed. AFR rate 25 also passed with a minor
descriptive upper-tail warning: subset q99.5 was 0.05797 versus a 0.05683
reference limit, while duplicate/MZ and first-degree counts were zero and all
frequency, LD, eigenstructure, PCA, and round-trip gates passed. Because this
limit came from 20 subsets of 200 chromosome-22 samples, q99.5 is descriptive
and is not an isolated hard biological boundary.

## Production results

Both canonical filesets contain exactly 5,000 individuals and 50,000 markers.
All calls are nonmissing 0/1/2 hard calls, marker identity is exact, individual
IDs are unique, and each SNP-major BED is exactly 62,500,003 bytes.

| Metric | EUR | AFR |
|---|---:|---:|
| MAF correlation | 0.999685 | 0.999650 |
| absolute mean-MAF difference | 0.000141 | 0.000062 |
| absolute heterozygosity difference | 0.002554 | 0.002103 |
| mean-r2 ratio, 0-10 kb | 0.9864 | 0.9895 |
| mean-r2 ratio, 10-50 kb | 0.9819 | 0.9907 |
| mean-r2 ratio, 50-100 kb | 0.9874 | 1.0022 |
| median LD-score ratio | 1.0161 | 1.0221 |
| effective-rank ratio | 1.0047 | 1.0068 |
| largest-eigenvalue ratio | 0.9087 | 0.8402 |
| maximum PCA centroid RMS z | 0.0533 | 0.0401 |
| duplicate/MZ pairs | 0 | 0 |
| first-degree pairs | 0 | 0 |
| second-degree pairs, descriptive | 3,079 | 167 |

AFR's synthetic-subset q99.5 median was 0.05775 versus the matched-reference
97.5% bound 0.05474, so the upper-tail warning persisted. Its full-panel q99.5
was 0.04952, with zero MZ and first-degree pairs and no corroborating close-pair
signal. The warning was therefore retained but did not fail production.

Full-panel frequency and relatedness summaries were combined with five
deterministic synthetic 200-person subsets and 20 matched reference subsets on
a deterministic 10,000-marker validation selection. LD and eigenstructure did
not require a dense 50,000-marker correlation matrix. Staged PLINK files passed
exact genio and focused qgg read-back before the atomic dual installation.

The complete build and validation took 1,063.2 seconds on this machine.
Production peak memory was not instrumented; the preceding scaling run observed
approximately 840 MB working set. Software was R 4.4.1, genio 1.1.2, qgg 1.1.6,
and PLINK 2.0.0-a.7.3 (8 Aug 2026).

Historical comparison: sim1000G preserved MAF but failed local-LD fidelity;
HAPNEST was not evaluated. Neither is recommended for these production panels.

Reproduction command:

```text
Rscript scripts/build_human_mosaic.R --mode=full --ancestry=both
```

After an interrupted build, rerun with `--resume`. See
`human_validation_summary.csv`, `human_validation_acceptance.csv`,
`human_validation_relatedness.csv`, `human_validation_runtime.csv`,
`human_validation_donor_summary.csv`, and `human_validation_checksums.csv` for
compact machine-readable evidence.

# Human 1000 Genomes EUR LD panel

Status: `available_validated`.

`human_1000g_eur.{bed,bim,fam}` contains 5,000 synthetic individuals and
50,000 biallelic chromosome 22 SNPs in PLINK1 SNP-major format. Coordinates are
GRCh37. Marker IDs, alleles, positions, and order come from the final 20130502
1000 Genomes Phase 3 phased reference and are shared exactly with the AFR
panel. Donors are CEU, GBR, FIN, IBS, and TSI.

The reference-conditioned base-R haplotype mosaic uses copying rate 20 and seed
20960818. Synthetic haplotypes copy only same-ancestry phased donors, without
mutation, frequency adjustment, ancestry mixing, imputation, or individual
filtering. These are not original 1000 Genomes participants.

Production validation passed all hard gates. File sizes are BED 62,500,003
bytes, BIM 2,363,322 bytes, and FAM 215,000 bytes. SHA-256:

- BED: `39e8ff00f11ebafe8d39e3d602b08df21dfd38f5972af41176d39261bdef6dce`
- BIM: `4d1d305a3761859ff3623498eaec0c3e6411a6c64f750b7b1e0a66f022127ca9`
- FAM: `2c23d281118649e85befdcb747f80816525efc0dee7dc0ec9400fc5d45678f48`

The original phased VCF SHA-256 is
`a90c16c4ff2b3196476d506ae13cb3047fae8670163c7c932c4b0239aef3daf5`;
the population-panel SHA-256 is
`b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858`;
and the deterministic marker-list SHA-256 is
`d8cc5e061c3ed2a6b9bbf06eb057833b9f96eabf846f7c1d7f53431b24b217e6`.
Software was R 4.4.1, genio 1.1.2, qgg 1.1.6, and PLINK
2.0.0-a.7.3 (8 Aug 2026).

Reproduce both panels atomically with:

```text
Rscript scripts/build_human_mosaic.R --mode=full --ancestry=both
```

See `validation/human_validation_report.md` and `datasets.csv` for provenance
and validation metrics.

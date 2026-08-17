# Protected independent human fixture

This existing fileset is referred to by downstream projects as
`human_independent`, but the repository directory is named
`simulated_human_data` and the prefix is `human`. The discrepancy is documented
rather than renamed. Do not move, rename, regenerate, or modify these files.

The PLINK files contain exactly 5,000 individuals × 50,000 markers, 0/1/2 hard
calls, and no missing genotypes. There are 5,000 unique individual IDs and
50,000 unique marker IDs across chromosomes 1–22. BIM files contain physical
and genetic-map positions and nucleotide alleles. FAM parents and sex are
unknown (`0`) and its phenotype field is `0`.

Auxiliary files have 5,000 × 14 (`human.covar`), 5,000 × 3 (`human.pheno`),
and 50,000 × 2 (`human.info`) fields. Their biological definitions, source,
generation method, genome build, and seed are not recorded in this repository;
this README does not invent them.

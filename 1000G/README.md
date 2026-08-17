# Existing small 1000 Genomes fixtures

These files are retained unchanged.

- `sample_chr1` is 489 individuals × 1,000 markers, but its BIM chromosome is
  21 (not chromosome 1). Positions span 9,412,099–11,043,740 bp.
- `sample_chr2` is 489 individuals × 1,000 markers, but its BIM chromosome is
  22 (not chromosome 2). Positions span 16,050,840–16,860,036 bp.

Both filesets contain exact 0/1/2 hard calls with no missing genotypes. They
have unique rsIDs, physical positions, genetic-map values, and nucleotide
alleles. Their identical FAM files contain 489 unique 1000 Genomes sample IDs;
FID equals IID, parents and sex are unknown (`0`), and phenotype is `-9`.

The extraction source, release, population selection, and meaning of the
legacy prefix names are not recorded in this repository, so they are not
inferred here. Some chromosome 21 genetic-map values are negative; they are
reported as stored and are not corrected without provenance.

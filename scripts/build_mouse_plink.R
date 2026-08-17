#!/usr/bin/env Rscript

# Build the canonical mouse PLINK1 fileset from the existing RDS objects.
# Unknown marker coordinates and alleles are represented by explicit PLINK
# placeholders; see mouse_data/README.md.

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

required_packages <- c("genio", "qgg", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "qgdata.Rproj")) &&
        dir.exists(file.path(path, "mouse_data"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the qgdata repository root")
    path <- parent
  }
}

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

root <- find_repo_root()
mouse_dir <- file.path(root, "mouse_data")
out_prefix <- file.path(mouse_dir, "mouse")
outputs <- paste0(out_prefix, c(".bed", ".bim", ".fam"))

if (any(file.exists(outputs)) && !force) {
  stop("Mouse PLINK output already exists; use --force to rebuild: ",
       paste(outputs[file.exists(outputs)], collapse = ", "))
}

raw <- as.matrix(readRDS(file.path(mouse_dir, "genotypes.rds")))
imputed <- as.matrix(readRDS(file.path(mouse_dir, "genotypes_imputed.rds")))
ped <- readRDS(file.path(mouse_dir, "pedigree.rds"))

assert(is.numeric(raw), "genotypes.rds must be numeric")
assert(is.numeric(imputed), "genotypes_imputed.rds must be numeric")
assert(identical(dim(raw), dim(imputed)), "Raw and imputed dimensions differ")
assert(identical(dimnames(raw), dimnames(imputed)), "Raw and imputed IDs/order differ")
assert(all(is.na(raw) | raw %in% 0:2), "genotypes.rds contains values outside 0/1/2/NA")
assert(all(is.finite(imputed)), "genotypes_imputed.rds contains missing/non-finite values")
assert(all(imputed %in% 0:2),
       "genotypes_imputed.rds is not exact 0/1/2 hard calls; refusing to round dosages")
assert(identical(raw[!is.na(raw)], imputed[!is.na(raw)]),
       "Imputed calls alter observed genotypes")
assert(all(c("id", "sire", "dam", "family", "sex") %in% names(ped)),
       "pedigree.rds lacks required pedigree columns")
assert(!anyDuplicated(ped$id), "pedigree.rds has duplicated IDs")
assert(identical(as.character(ped$id), rownames(imputed)),
       "Pedigree IDs do not exactly match genotype row IDs/order")

# genio expects markers x individuals and counts of the BIM alternative allele.
X <- t(imputed)
rownames(X) <- colnames(imputed)
colnames(X) <- rownames(imputed)

# No chromosome, coordinate, map, or nucleotide allele metadata exists in any
# mouse RDS object. Zero denotes unknown PLINK fields; 1/2 are non-biological
# allele labels matching the source dosage coding.
bim <- data.frame(
  chr = rep(0L, nrow(X)),
  id = rownames(X),
  posg = rep(0, nrow(X)),
  pos = rep(0L, nrow(X)),
  ref = rep("1", nrow(X)),
  alt = rep("2", nrow(X)),
  stringsAsFactors = FALSE
)

sex_code <- unname(c(Male = 1L, Female = 2L)[as.character(ped$sex)])
assert(!anyNA(sex_code), "Unexpected sex value in pedigree.rds")
fam <- data.frame(
  fam = as.character(ped$family),
  id = as.character(ped$id),
  pat = as.character(ped$sire),
  mat = as.character(ped$dam),
  sex = sex_code,
  pheno = rep(-9, nrow(ped)),
  stringsAsFactors = FALSE
)

tmp_dir <- tempfile("qgdata-mouse-plink-")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
tmp_prefix <- file.path(tmp_dir, "mouse")

genio::write_plink(tmp_prefix, X, bim = bim, fam = fam, verbose = FALSE)

# Independent writer read-back before qgg validation.
back_genio <- genio::read_plink(tmp_prefix, verbose = FALSE)
assert(identical(dim(back_genio$X), dim(X)), "genio read-back dimensions differ")
assert(identical(rownames(back_genio$X), rownames(X)), "genio marker order differs")
assert(identical(colnames(back_genio$X), colnames(X)), "genio individual order differs")
assert(identical(is.na(back_genio$X), is.na(X)), "genio missing-value pattern differs")
assert(isTRUE(all.equal(unname(back_genio$X), unname(X),
                        check.attributes = FALSE, tolerance = 0)),
       "genio genotype read-back differs")

# qgg uses chr as a fileset index, not the BIM chromosome label.
qgg_tmp <- file.path(tmp_dir, "qgg")
dir.create(qgg_tmp)
old_wd <- setwd(qgg_tmp)
on.exit(setwd(old_wd), add = TRUE)
glist <- qgg::gprep(
  study = "qgdata_mouse_roundtrip",
  bedfiles = paste0(tmp_prefix, ".bed"),
  bimfiles = paste0(tmp_prefix, ".bim"),
  famfiles = paste0(tmp_prefix, ".fam"),
  overwrite = TRUE
)
back_qgg <- qgg::getG(Glist = glist, chr = 1, impute = FALSE, scale = FALSE)
setwd(old_wd)

assert(identical(dim(back_qgg), dim(imputed)), "qgg read-back dimensions differ")
assert(identical(rownames(back_qgg), rownames(imputed)), "qgg individual order differs")
assert(identical(colnames(back_qgg), colnames(imputed)), "qgg marker order differs")
assert(identical(is.na(back_qgg), is.na(imputed)), "qgg missing-value pattern differs")
assert(isTRUE(all.equal(unname(back_qgg), unname(imputed),
                        check.attributes = FALSE, tolerance = 0)),
       "qgg genotype read-back differs")

fam_back <- utils::read.table(paste0(tmp_prefix, ".fam"), header = FALSE,
                              colClasses = "character", stringsAsFactors = FALSE)
assert(identical(fam_back[[1]], fam$fam), "FAM family IDs differ")
assert(identical(fam_back[[2]], fam$id), "FAM individual IDs differ")
assert(identical(fam_back[[3]], fam$pat), "FAM paternal IDs differ")
assert(identical(fam_back[[4]], fam$mat), "FAM maternal IDs differ")
assert(identical(as.integer(fam_back[[5]]), fam$sex), "FAM sex differs")
assert(all(fam_back[[6]] == "-9"), "FAM phenotype is not -9")

# Copy only after every validation passes. Existing outputs are overwritten only
# under the explicit --force flag checked above.
for (ext in c(".bed", ".bim", ".fam")) {
  src <- paste0(tmp_prefix, ext)
  dst <- paste0(out_prefix, ext)
  ok <- file.copy(src, dst, overwrite = force)
  assert(ok, paste("Failed to install", dst))
}

sha256 <- vapply(outputs, digest::digest, character(1), algo = "sha256", file = TRUE)
report <- c(
  "dataset=mouse",
  "source=mouse_data/genotypes_imputed.rds",
  paste0("individuals=", nrow(imputed)),
  paste0("markers=", ncol(imputed)),
  paste0("raw_missing_calls=", sum(is.na(raw))),
  "imputed_missing_calls=0",
  "imputed_values=exact_integer_0_1_2",
  "observed_raw_calls_unchanged=TRUE",
  "genio_roundtrip_exact=TRUE",
  "qgg_gprep_getG_roundtrip_exact=TRUE",
  "pedigree_roundtrip_exact=TRUE",
  paste0(basename(outputs), "_sha256=", sha256),
  paste0("R=", R.version.string),
  paste0("genio=", as.character(utils::packageVersion("genio"))),
  paste0("qgg=", as.character(utils::packageVersion("qgg")))
)
writeLines(report, file.path(mouse_dir, "mouse_validation.txt"))
capture.output(sessionInfo(), file = file.path(mouse_dir, "mouse_sessionInfo.txt"))
cat(paste(report, collapse = "\n"), "\n")

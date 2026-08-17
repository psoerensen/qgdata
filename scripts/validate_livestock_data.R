#!/usr/bin/env Rscript

# Focused validation for the canonical cattle and population-specific pig
# datasets. Large matrices are never materialized; full PED equivalence is
# checked one individual at a time.

required <- c("data.table", "digest", "qgg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "qgdata.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find qgdata repository root")
    path <- parent
  }
}
sha256 <- function(path) unname(digest::digest(path, algo = "sha256", file = TRUE))
root <- find_root()
cache <- file.path(root, "cache", "livestock")
plink2 <- file.path(root, "cache", "tools", "plink2_20260808", "plink2.exe")
assert(file.exists(plink2), paste("Missing cached PLINK2:", plink2))
work_root <- file.path(cache, paste0("validation-", Sys.getpid()))
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_root, recursive = TRUE, force = TRUE), add = TRUE)

run_plink2 <- function(arguments, label) {
  output <- system2(plink2, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L)
    stop("PLINK2 failed during ", label, ":\n", paste(output, collapse = "\n"))
  invisible(output)
}
plink_version <- paste(run_plink2("--version", "version check"), collapse = " ")

read_fam <- function(prefix) {
  x <- utils::read.table(paste0(prefix, ".fam"), header = FALSE,
                         colClasses = "character", stringsAsFactors = FALSE)
  assert(ncol(x) == 6L, paste("Invalid FAM:", prefix))
  names(x) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
  x
}
read_bim <- function(prefix) {
  x <- data.table::fread(paste0(prefix, ".bim"), header = FALSE,
                         colClasses = "character", data.table = FALSE)
  assert(ncol(x) == 6L, paste("Invalid BIM:", prefix))
  names(x) <- c("CHR", "ID", "CM", "POS", "A1", "A2")
  x
}
expected_bed_bytes <- function(n, m) 3 + as.double(m) * ceiling(n / 4)

validate_plink_structure <- function(prefix, expected_n, expected_m, label,
                                     plink_options = character()) {
  files <- paste0(prefix, c(".bed", ".bim", ".fam"))
  assert(all(file.exists(files)), paste("Incomplete PLINK fileset:", label))
  con <- file(paste0(prefix, ".bed"), "rb")
  magic <- readBin(con, "raw", 3L); close(con)
  assert(identical(as.integer(magic), c(108L, 27L, 1L)),
         paste(label, "BED is not SNP-major PLINK1"))
  fam <- read_fam(prefix); bim <- read_bim(prefix)
  assert(nrow(fam) == expected_n && nrow(bim) == expected_m,
         paste(label, "dimensions differ"))
  assert(!anyDuplicated(paste(fam$FID, fam$IID, sep = "\r")),
         paste(label, "has duplicate FID/IID pairs"))
  assert(!anyDuplicated(bim$ID), paste(label, "has duplicate marker IDs"))
  assert(all(nzchar(bim$CHR)) && all(is.finite(as.numeric(bim$POS))) &&
           all(as.numeric(bim$POS) >= 0) && all(nzchar(bim$A1)) && all(nzchar(bim$A2)),
         paste(label, "has invalid BIM fields"))
  assert(all(fam$SEX %in% c("0", "1", "2")) && all(nzchar(fam$PAT)) &&
           all(nzchar(fam$MAT)), paste(label, "has invalid FAM pedigree/sex fields"))
  actual_size <- unname(file.info(paste0(prefix, ".bed"))$size)
  assert(actual_size == expected_bed_bytes(expected_n, expected_m),
         paste(label, "BED byte size differs"))

  miss_prefix <- file.path(work_root, paste0(gsub("[^A-Za-z0-9]", "_", label), "_missing"))
  run_plink2(c(plink_options, "--bfile", prefix, "--missing", "--out", miss_prefix),
             paste(label, "missingness"))
  smiss <- data.table::fread(paste0(miss_prefix, ".smiss"), data.table = FALSE)
  names(smiss) <- sub("^#", "", names(smiss))
  missing_name <- intersect(c("MISSING_CT", "N_MISS"), names(smiss))
  observed_name <- intersect(c("OBS_CT", "N_GENO"), names(smiss))
  assert(length(missing_name) == 1L && length(observed_name) == 1L,
         paste(label, "unexpected PLINK missingness schema"))
  missing_calls <- sum(smiss[[missing_name]])
  observed_calls <- sum(smiss[[observed_name]])
  # PLINK's OBS_CT is sex-aware on X/Y; it can therefore be smaller than the
  # rectangular BED cell count when female samples and Y markers are present.
  assert(observed_calls > 0 && observed_calls <= expected_n * expected_m,
         paste(label, "PLINK reported an invalid observed call count"))
  list(fam = fam, bim = bim, missing_calls = missing_calls,
       plink_options = plink_options,
       missing_fraction = missing_calls / observed_calls, bed_bytes = actual_size)
}

focused_roundtrip <- function(prefix, structure, label) {
  focus_dir <- file.path(work_root, paste0("focus_", gsub("[^A-Za-z0-9]", "_", label)))
  dir.create(focus_dir)
  sample_index <- unique(round(seq(1, nrow(structure$fam), length.out = min(50L, nrow(structure$fam)))))
  marker_index <- unique(round(seq(1, nrow(structure$bim), length.out = min(100L, nrow(structure$bim)))))
  keep <- file.path(focus_dir, "keep.txt")
  extract <- file.path(focus_dir, "extract.txt")
  utils::write.table(structure$fam[sample_index, c("FID", "IID")], keep,
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(structure$bim$ID[marker_index], extract)
  out <- file.path(focus_dir, "focused")
  run_plink2(c(structure$plink_options, "--bfile", prefix,
               "--keep", keep, "--extract", extract,
               "--make-bed", "--out", out), paste(label, "focused extraction"))

  qgg_dir <- file.path(focus_dir, "qgg"); dir.create(qgg_dir)
  old <- setwd(qgg_dir); on.exit(setwd(old), add = TRUE)
  glist <- qgg::gprep(study = paste0("livestock_", gsub("[^A-Za-z0-9]", "_", label)),
                      bedfiles = paste0(out, ".bed"), bimfiles = paste0(out, ".bim"),
                      famfiles = paste0(out, ".fam"), overwrite = TRUE)
  setwd(old)
  qgg_back <- qgg::getG(Glist = glist, chr = 1, impute = FALSE, scale = FALSE)
  assert(all(is.na(qgg_back) | qgg_back %in% 0:2),
         paste(label, "qgg returned invalid genotype values"))

  genio_status <- "not_installed"
  if (requireNamespace("genio", quietly = TRUE)) {
    genio_back <- t(genio::read_plink(out, verbose = FALSE)$X)
    assert(identical(dim(genio_back), dim(qgg_back)) &&
             identical(rownames(genio_back), rownames(qgg_back)) &&
             identical(colnames(genio_back), colnames(qgg_back)) &&
             identical(is.na(genio_back), is.na(qgg_back)) &&
             isTRUE(all.equal(unname(genio_back), unname(qgg_back), tolerance = 0,
                              check.attributes = FALSE)),
           paste(label, "focused genio/qgg read-backs differ"))
    genio_status <- "exact"
  }
  list(qgg = "valid_0_1_2_NA", genio = genio_status,
       individuals = nrow(qgg_back), markers = ncol(qgg_back))
}

compare_ped_streaming <- function(source_ped, prefix, expected_markers, label) {
  export_prefix <- file.path(work_root, paste0(gsub("[^A-Za-z0-9]", "_", label), "_export"))
  run_plink2(c("--bfile", prefix, "--export", "ped", "--out", export_prefix),
             paste(label, "PED round-trip export"))
  source_con <- file(source_ped, "rt"); export_con <- file(paste0(export_prefix, ".ped"), "rt")
  on.exit({ close(source_con); close(export_con) }, add = TRUE)
  row <- 0L
  repeat {
    source_line <- readLines(source_con, n = 1L, warn = FALSE)
    export_line <- readLines(export_con, n = 1L, warn = FALSE)
    if (!length(source_line) || !length(export_line)) {
      assert(!length(source_line) && !length(export_line),
             paste(label, "PED row counts differ after round-trip"))
      break
    }
    row <- row + 1L
    source_tokens <- strsplit(trimws(source_line), "[[:space:]]+")[[1L]]
    export_tokens <- strsplit(trimws(export_line), "[[:space:]]+")[[1L]]
    assert(length(source_tokens) == 6L + 2L * expected_markers &&
             length(export_tokens) == length(source_tokens),
           paste(label, "PED field count differs at row", row))
    assert(identical(source_tokens[1:6], export_tokens[1:6]),
           paste(label, "PED metadata differs at row", row))
    source_genotype <- as.integer(source_tokens[-(1:6)])
    export_genotype <- as.integer(export_tokens[-(1:6)])
    assert(!anyNA(source_genotype) && !anyNA(export_genotype) &&
             all(source_genotype %in% 0:2) && all(export_genotype %in% 0:2),
           paste(label, "PED has invalid allele codes at row", row))
    odd <- seq.int(1L, length(source_genotype), 2L); even <- odd + 1L
    source_missing <- source_genotype[odd] == 0L | source_genotype[even] == 0L
    export_missing <- export_genotype[odd] == 0L | export_genotype[even] == 0L
    assert(identical(source_missing, export_missing) &&
             all(source_genotype[odd] + source_genotype[even] ==
                   export_genotype[odd] + export_genotype[even]) &&
             all(source_genotype[odd] * source_genotype[even] ==
                   export_genotype[odd] * export_genotype[even]),
           paste(label, "genotype differs after allele-orientation-aware round-trip at row", row))
  }
  row
}

trait_statistics <- function(x, dataset, population = "all") {
  id_columns <- if (dataset == "cattle") c("FID", "IID") else
    c("FID", "IID", "population")
  traits <- setdiff(names(x), id_columns)
  do.call(rbind, lapply(traits, function(trait) {
    values <- suppressWarnings(as.numeric(x[[trait]]))
    nonmissing <- values[!is.na(values)]
    data.frame(dataset = dataset, population = population, trait = trait,
      total = length(values), nonmissing = length(nonmissing), missing = sum(is.na(values)),
      mean = if (length(nonmissing)) mean(nonmissing) else NA_real_,
      sd = if (length(nonmissing) > 1L) stats::sd(nonmissing) else NA_real_,
      min = if (length(nonmissing)) min(nonmissing) else NA_real_,
      max = if (length(nonmissing)) max(nonmissing) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

write_validation <- function(directory, lines) {
  path <- file.path(directory, "validation.txt")
  temporary <- paste0(path, ".partial-", Sys.getpid())
  writeLines(lines, temporary)
  assert(file.rename(temporary, path), paste("Could not write", path))
}

# Cattle validation.
cattle_source <- file.path(cache, "source", "cattle")
cattle_dir <- file.path(root, "cattle_milk_lipids")
cattle_prefix <- file.path(cattle_dir, "cattle_milk_lipids")
assert(dir.exists(cattle_source) && dir.exists(cattle_dir), "Cattle source/output is missing")
cattle_source_fam <- read_fam(file.path(cattle_source, "GenotypesHD"))
cattle_source_bim <- read_bim(file.path(cattle_source, "GenotypesHD"))
cattle <- validate_plink_structure(cattle_prefix, nrow(cattle_source_fam),
                                   nrow(cattle_source_bim), "cattle", "--cow")
for (extension in c("bed", "bim", "fam"))
  assert(identical(sha256(file.path(cattle_source, paste0("GenotypesHD.", extension))),
                   sha256(paste0(cattle_prefix, ".", extension))),
         paste("Cattle installed", extension, "is not byte-identical to source"))

cattle_pheno_path <- file.path(cattle_dir, "phenotypes.csv")
assert(identical(sha256(file.path(cattle_source, "PhenoPLs.csv")), sha256(cattle_pheno_path)),
       "Cattle phenotype CSV is not byte-identical to source")
cattle_pheno <- data.table::fread(cattle_pheno_path, na.strings = "NA",
                                  colClasses = "character", data.table = FALSE)
cattle_pheno_ids <- paste(cattle_pheno$FID, cattle_pheno$IID, sep = "\r")
cattle_fam_ids <- paste(cattle$fam$FID, cattle$fam$IID, sep = "\r")
assert(nrow(cattle_pheno) == nrow(cattle$fam) && ncol(cattle_pheno) == 58L &&
         !anyDuplicated(cattle_pheno_ids) && setequal(cattle_pheno_ids, cattle_fam_ids),
       "Cattle phenotype IDs do not match FAM")
source_fixed <- data.table::fread(file.path(cattle_source, "FixedEffectPLs.txt"),
                                  colClasses = "character", data.table = FALSE)
installed_fixed <- data.table::fread(file.path(cattle_dir, "fixed_effects.csv"),
                                     colClasses = "character", data.table = FALSE)
fixed_ids <- paste(installed_fixed$FID, installed_fixed$IID, sep = "\r")
assert(identical(source_fixed, installed_fixed) && !anyDuplicated(fixed_ids) &&
         setequal(fixed_ids, cattle_fam_ids),
       "Cattle fixed effects changed or IDs differ")
cattle_focus <- focused_roundtrip(cattle_prefix, cattle, "cattle")
cattle_stats <- trait_statistics(cattle_pheno, "cattle")

# Pig validation.
pig_source <- file.path(cache, "source", "pig", "Primary_data", "GWAS")
pig_dir <- file.path(root, "pig_blood_lipids")
assert(dir.exists(pig_source) && dir.exists(pig_dir), "Pig source/output is missing")
definitions <- list(
  laiwu = list(source = "Laiwu", population = "Laiwu"),
  erhualian = list(source = "EHL", population = "Erhualian"),
  dly = list(source = "DLY", population = "Duroc_x_Landrace_x_Yorkshire")
)
pig_results <- list(); pig_focus <- list(); pig_roundtrip_rows <- integer()
for (name in names(definitions)) {
  definition <- definitions[[name]]
  source_map_path <- file.path(pig_source, paste0(definition$source, "_60K_gens.map"))
  source_ped_path <- file.path(pig_source, paste0(definition$source, "_60K_gens.ped"))
    # read.table() is intentional: one archived marker ID begins with '#',
    # which fread treats as a comment and would silently omit.
    source_map <- utils::read.table(source_map_path, header = FALSE,
                                    colClasses = "character", quote = "",
                                    comment.char = "", stringsAsFactors = FALSE)
  source_ids <- data.table::fread(source_ped_path, header = FALSE, select = 1:6,
                                  colClasses = "character", data.table = FALSE,
                                  showProgress = FALSE)
  prefix <- file.path(pig_dir, name, name)
  result <- validate_plink_structure(prefix, nrow(source_ids), nrow(source_map),
                                     paste("pig", name))
    expected_map <- source_map[, 1:4, drop = FALSE]
    # The archived files use PLINK's documented numeric aliases 23/24 for
    # X/Y.  PLINK2 emits the canonical X/Y labels in BIM while retaining
    # every marker, coordinate, and its order.
    expected_map[[1L]][expected_map[[1L]] == "23"] <- "X"
    expected_map[[1L]][expected_map[[1L]] == "24"] <- "Y"
    assert(identical(unname(expected_map),
                     unname(result$bim[, c("CHR", "ID", "CM", "POS")])),
           paste("Pig", name, "MAP fields/order changed"))
  assert(identical(unname(source_ids), unname(result$fam)),
         paste("Pig", name, "FAM metadata/order changed"))
    assert(all(result$bim$A1 %in% c(".", "1", "2")) &&
             all(result$bim$A2 %in% c(".", "1", "2")),
           paste("Pig", name, "has alleles not present in recoded source"))
  pig_roundtrip_rows[[name]] <- compare_ped_streaming(source_ped_path, prefix,
                                                      nrow(source_map), paste("pig", name))
  pig_focus[[name]] <- focused_roundtrip(prefix, result, paste("pig", name))
  pig_results[[name]] <- result
}

pig_pheno <- data.table::fread(file.path(pig_dir, "phenotypes.csv"),
                               colClasses = "character", na.strings = NULL,
                               data.table = FALSE, check.names = FALSE)
traits <- c("TCHOL", "TG", "HDL-C", "LDL-C", "HDL-C/LDL-C", "AI")
assert(identical(names(pig_pheno), c("FID", "IID", "population", traits)),
       "Canonical pig phenotype columns differ")
expected_pig_ids <- do.call(rbind, lapply(names(definitions), function(name) {
  x <- pig_results[[name]]$fam
  data.frame(FID = x$FID, IID = x$IID, population = definitions[[name]]$population,
             stringsAsFactors = FALSE)
}))
assert(identical(pig_pheno[, c("FID", "IID", "population")], expected_pig_ids),
       "Pig phenotype IDs/order do not exactly match population FAM files")

phenotype_only <- integer(); genotype_only <- integer()
for (name in names(definitions)) {
  definition <- definitions[[name]]
  path <- file.path(pig_source, paste0(definition$source, "_60K_phens.txt"))
  source_pheno <- utils::read.delim(path, header = TRUE, colClasses = "character",
    check.names = FALSE, quote = "", comment.char = "", na.strings = NULL,
    stringsAsFactors = FALSE)
  names(source_pheno)[names(source_pheno) == "HDL.C"] <- "HDL-C"
  names(source_pheno)[names(source_pheno) == "LDL.C"] <- "LDL-C"
  canonical <- pig_pheno[pig_pheno$population == definition$population, , drop = FALSE]
  index <- match(canonical$IID, source_pheno$id)
  expected <- matrix("NA", nrow(canonical), length(traits), dimnames = list(NULL, traits))
  matched <- !is.na(index)
  expected[matched, ] <- as.matrix(source_pheno[index[matched], traits, drop = FALSE])
  assert(identical(unname(as.matrix(canonical[, traits, drop = FALSE])), unname(expected)),
         paste("Pig", name, "phenotype values changed"))
  phenotype_only[[name]] <- sum(!source_pheno$id %in% canonical$IID)
  genotype_only[[name]] <- sum(!canonical$IID %in% source_pheno$id)
}
pig_stats <- do.call(rbind, lapply(unique(pig_pheno$population), function(population)
  trait_statistics(pig_pheno[pig_pheno$population == population, , drop = FALSE],
                   "pig", population)))

format_stats <- function(x) apply(x, 1L, function(row) paste0(
  "trait=", row[["trait"]], ";population=", row[["population"]],
  ";total=", row[["total"]], ";nonmissing=", row[["nonmissing"]],
  ";missing=", row[["missing"]], ";mean=", row[["mean"]],
  ";sd=", row[["sd"]], ";min=", row[["min"]], ";max=", row[["max"]]))

data_file_lines <- function(directory, exclude = "validation.txt") {
  files <- list.files(directory, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE & basename(files) != exclude]
  relative <- substring(normalizePath(files, winslash = "/"),
                        nchar(normalizePath(directory, winslash = "/")) + 2L)
  info <- file.info(files)
  assert(all(info$size < 100000000), paste("File at/above 100 MB in", directory))
  paste0("file=", relative, ";bytes=", info$size,
         ";sha256=", vapply(files, sha256, character(1)))
}

archive_cattle <- file.path(cache, "outer_extracted", "cattle",
                            "Polar_Lipids_Phenotypes_HDGenotypes.zip")
archive_pig <- file.path(cache, "outer_extracted", "pig", "Primary_data.zip")
common_versions <- c(paste0("R=", R.version.string), paste0("PLINK2=", plink_version),
                     paste0("qgg=", as.character(utils::packageVersion("qgg"))),
                     paste0("genio=", if (requireNamespace("genio", quietly = TRUE))
                       as.character(utils::packageVersion("genio")) else "not_installed"))

cattle_lines <- c(
  "status=PASS", "dataset=cattle_milk_lipids", "license=CC0-1.0",
  "doi=https://doi.org/10.5061/dryad.bcc2fqzph",
  "version_url=https://datadryad.org/api/v2/versions/368250/download",
  paste0("source_archive_bytes=", file.info(archive_cattle)$size),
  paste0("source_archive_sha256=", sha256(archive_cattle)),
  paste0("individuals=", nrow(cattle$fam)), paste0("markers=", nrow(cattle$bim)),
  paste0("missing_genotypes=", cattle$missing_calls),
  paste0("missing_genotype_fraction=", cattle$missing_fraction),
  "source_plink_byte_identical=TRUE", "phenotypes_byte_identical=TRUE",
  "fixed_effect_values_exact=TRUE", "phenotype_ids_match=TRUE",
  paste0("qgg_roundtrip=", cattle_focus$qgg),
  paste0("genio_roundtrip=", cattle_focus$genio), common_versions,
  format_stats(cattle_stats), data_file_lines(cattle_dir))
write_validation(cattle_dir, cattle_lines)

pig_lines <- c(
  "status=PASS", "dataset=pig_blood_lipids", "license=CC0-1.0",
  "doi=https://doi.org/10.5061/dryad.4gh70",
  "version_url=https://datadryad.org/api/v2/versions/25357/download",
  paste0("source_archive_bytes=", file.info(archive_pig)$size),
  paste0("source_archive_md5=", unname(tools::md5sum(archive_pig))),
  paste0("source_archive_sha256=", sha256(archive_pig)),
  paste0("population=", names(pig_results), ";individuals=",
         vapply(pig_results, function(x) nrow(x$fam), integer(1)), ";markers=",
         vapply(pig_results, function(x) nrow(x$bim), integer(1)),
         ";missing_genotypes=", vapply(pig_results, function(x) x$missing_calls, numeric(1)),
         ";missing_fraction=", vapply(pig_results, function(x) x$missing_fraction, numeric(1))),
  paste0("population=", names(pig_results), ";ped_roundtrip_rows=", pig_roundtrip_rows,
         ";genotype_equivalence=TRUE;qgg_roundtrip=",
         vapply(pig_focus, function(x) x$qgg, character(1)), ";genio_roundtrip=",
         vapply(pig_focus, function(x) x$genio, character(1))),
  paste0("population=", names(definitions), ";phenotype_only_source_rows=", phenotype_only,
         ";genotype_only_rows=", genotype_only),
  "phenotype_ids_match=TRUE", "phenotype_values_exact=TRUE",
  "populations_deliberately_separate=TRUE", common_versions,
  format_stats(pig_stats), data_file_lines(pig_dir))
write_validation(pig_dir, pig_lines)

# Validation files themselves are checked for the GitHub size guard after they
# are finalized; a file cannot contain its own stable checksum.
all_final <- c(list.files(cattle_dir, recursive = TRUE, full.names = TRUE),
               list.files(pig_dir, recursive = TRUE, full.names = TRUE))
assert(all(file.info(all_final)$size < 100000000), "A final livestock file is at/above 100 MB")

cat("status=PASS\n")
cat("cattle=", nrow(cattle$fam), "x", nrow(cattle$bim), "\n", sep = "")
for (name in names(pig_results))
  cat("pig_", name, "=", nrow(pig_results[[name]]$fam), "x",
      nrow(pig_results[[name]]$bim), "\n", sep = "")
cat("cattle_traits=56\n")
cat("pig_traits=6\n")
cat("pig_ped_roundtrip_equivalence=TRUE\n")
cat("focused_qgg_roundtrip=TRUE\n")
cat("focused_genio_roundtrip=",
    if (requireNamespace("genio", quietly = TRUE)) "TRUE" else "not_installed", "\n", sep = "")
cat("all_final_files_below_100MB=TRUE\n")

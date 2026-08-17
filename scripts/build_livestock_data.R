#!/usr/bin/env Rscript

# Reproducible livestock dataset builder. Source archives and all PED/MAP
# intermediates remain under ignored cache/.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[length(hit)]]) else default
}
species <- match.arg(tolower(arg_value("species", "both")),
                     c("cattle", "pig", "both"))
force <- "--force" %in% args

required <- c("data.table", "digest")
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
archive_dir <- file.path(cache, "archives")
outer_dir <- file.path(cache, "outer_extracted")
source_dir <- file.path(cache, "source")
plink2 <- file.path(root, "cache", "tools", "plink2_20260808", "plink2.exe")
assert(file.exists(plink2), paste("Missing cached PLINK2:", plink2))
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

sources <- list(
  cattle = list(
    version = "368250",
    url = "https://datadryad.org/api/v2/versions/368250/download",
    outer = "dryad_version_368250.zip",
    inner = "Polar_Lipids_Phenotypes_HDGenotypes.zip",
    inner_bytes = 16745878,
    inner_digest = "d92af1f1ddd63b401b18e88d6b65ed20cdc5c366cdff24fde16224c1d34f9590",
    digest_type = "sha256"
  ),
  pig = list(
    version = "25357",
    url = "https://datadryad.org/api/v2/versions/25357/download",
    outer = "dryad_version_25357.zip",
    inner = "Primary_data.zip",
    inner_bytes = 22332040,
    inner_digest = "51bcd5b1011b55379a0425bf5c0e2a2c",
    digest_type = "md5"
  )
)

download_with_retries <- function(url, destination) {
  if (file.exists(destination)) return(invisible(destination))
  options(timeout = 600)
  errors <- character()
  for (attempt in seq_len(3L)) {
    temporary <- paste0(destination, ".partial")
    if (file.exists(temporary)) unlink(temporary)
    result <- tryCatch({
      status <- suppressWarnings(utils::download.file(
        url, temporary, method = "libcurl", mode = "wb", quiet = FALSE,
        headers = c("User-Agent" = "qgdata/1.0", "X-API-Version" = "2.1.0")
      ))
      identical(status, 0L)
    }, error = function(e) {
      errors <<- c(errors, conditionMessage(e)); FALSE
    })
    if (isTRUE(result) && file.exists(temporary) && file.info(temporary)$size > 0) {
      assert(file.rename(temporary, destination), paste("Could not finalize", destination))
      return(invisible(destination))
    }
    if (file.exists(temporary)) unlink(temporary)
    if (attempt < 3L) Sys.sleep(2)
  }
  stop("Download failed after three attempts: ", url, "\n", paste(errors, collapse = "\n"))
}

prepare_source <- function(kind) {
  info <- sources[[kind]]
  outer <- file.path(archive_dir, info$outer)
  download_with_retries(info$url, outer)
  con <- file(outer, "rb")
  magic <- readBin(con, "raw", 2L)
  close(con)
  assert(identical(as.integer(magic), c(80L, 75L)),
         paste("Dryad version response is not ZIP:", outer))

  extracted_outer <- file.path(outer_dir, kind)
  inner <- file.path(extracted_outer, info$inner)
  if (!file.exists(inner)) {
    dir.create(extracted_outer, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(outer, exdir = extracted_outer)
  }
  candidates <- list.files(extracted_outer, pattern = paste0("^", info$inner, "$"),
                           recursive = TRUE, full.names = TRUE)
  assert(length(candidates) == 1L, paste("Expected one inner archive:", info$inner))
  inner <- candidates[[1L]]
  assert(file.info(inner)$size == info$inner_bytes,
         paste("Inner archive size mismatch:", info$inner))
  actual_digest <- if (info$digest_type == "sha256") sha256(inner) else
    unname(tools::md5sum(inner))
  assert(identical(tolower(actual_digest), info$inner_digest),
         paste("Inner archive checksum mismatch:", info$inner))

  extracted_source <- file.path(source_dir, kind)
  if (!dir.exists(extracted_source)) {
    dir.create(extracted_source, recursive = TRUE)
    utils::unzip(inner, exdir = extracted_source)
  }
  list(path = extracted_source, outer = outer, inner = inner,
       outer_bytes = unname(file.info(outer)$size), outer_sha256 = sha256(outer),
       inner_bytes = info$inner_bytes, inner_digest = actual_digest,
       inner_digest_type = info$digest_type, url = info$url, version = info$version)
}

run_plink2 <- function(arguments, label) {
  output <- system2(plink2, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L)
    stop("PLINK2 failed during ", label, ":\n", paste(output, collapse = "\n"))
  invisible(output)
}

copy_exact <- function(source, destination) {
  assert(file.copy(source, destination, overwrite = FALSE),
         paste("Could not copy", source))
  assert(identical(sha256(source), sha256(destination)),
         paste("Copied file differs:", destination))
}

build_cattle <- function(source, staging_root) {
  src <- source$path
  stage <- file.path(staging_root, "cattle_milk_lipids")
  dir.create(stage, recursive = TRUE)
  for (extension in c("bed", "bim", "fam"))
    copy_exact(file.path(src, paste0("GenotypesHD.", extension)),
               file.path(stage, paste0("cattle_milk_lipids.", extension)))
  copy_exact(file.path(src, "PhenoPLs.csv"), file.path(stage, "phenotypes.csv"))

  fixed <- data.table::fread(file.path(src, "FixedEffectPLs.txt"),
                             colClasses = "character", data.table = FALSE)
  assert(identical(names(fixed), c("FID", "IID", "Year_Batch")) && nrow(fixed) == 336L,
         "Unexpected cattle fixed-effect table")
  utils::write.table(fixed, file.path(stage, "fixed_effects.csv"), sep = ",",
                     quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
  invisible(source)
}

read_source_phenotype <- function(path) {
  utils::read.delim(path, header = TRUE, colClasses = "character",
                    check.names = FALSE, quote = "", comment.char = "",
                    na.strings = NULL, stringsAsFactors = FALSE)
}

build_pig <- function(source, staging_root) {
  gwas <- file.path(source$path, "Primary_data", "GWAS")
  assert(dir.exists(gwas), paste("Missing pig GWAS directory:", gwas))
  stage <- file.path(staging_root, "pig_blood_lipids")
  dir.create(stage, recursive = TRUE)
  definitions <- list(
    laiwu = list(prefix = "Laiwu", label = "Laiwu"),
    erhualian = list(prefix = "EHL", label = "Erhualian"),
    dly = list(prefix = "DLY", label = "Duroc_x_Landrace_x_Yorkshire")
  )
  traits <- c("TCHOL", "TG", "HDL-C", "LDL-C", "HDL-C/LDL-C", "AI")
  phenotype_rows <- list()

  for (name in names(definitions)) {
    definition <- definitions[[name]]
    source_prefix <- file.path(gwas, paste0(definition$prefix, "_60K_gens"))
    out_dir <- file.path(stage, name)
    dir.create(out_dir)
    out_prefix <- file.path(out_dir, name)
    run_plink2(c("--pedmap", source_prefix, "--make-bed", "--out", out_prefix),
               paste(name, "PED/MAP conversion"))
    unlink(paste0(out_prefix, ".log"))

    fam <- utils::read.table(paste0(out_prefix, ".fam"), header = FALSE,
                             colClasses = "character", stringsAsFactors = FALSE)
    names(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
    pheno_path <- file.path(gwas, paste0(definition$prefix, "_60K_phens.txt"))
    pheno <- read_source_phenotype(pheno_path)
    names(pheno)[names(pheno) == "HDL.C"] <- "HDL-C"
    names(pheno)[names(pheno) == "LDL.C"] <- "LDL-C"
    assert(all(c("id", traits) %in% names(pheno)),
           paste("Pig phenotype columns differ:", pheno_path))
    match_index <- match(fam$IID, pheno$id)
    values <- matrix("NA", nrow(fam), length(traits),
                     dimnames = list(NULL, traits))
    matched <- !is.na(match_index)
    values[matched, ] <- as.matrix(pheno[match_index[matched], traits, drop = FALSE])
    phenotype_rows[[name]] <- data.frame(
      FID = fam$FID, IID = fam$IID, population = definition$label,
      values, check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  combined <- do.call(rbind, phenotype_rows)
  rownames(combined) <- NULL
  assert(!anyDuplicated(paste(combined$population, combined$FID, combined$IID, sep = "\r")),
         "Combined pig phenotype table has duplicate genotype IDs")
  utils::write.table(combined, file.path(stage, "phenotypes.csv"), sep = ",",
                     quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
  invisible(source)
}

requested <- if (species == "both") c("cattle", "pig") else species
targets <- c(cattle = file.path(root, "cattle_milk_lipids"),
             pig = file.path(root, "pig_blood_lipids"))
existing <- requested[vapply(targets[requested], dir.exists, logical(1))]
if (length(existing) && !force)
  stop("Output already exists for ", paste(existing, collapse = ", "),
       "; use --force to replace it intentionally")

staging_root <- file.path(cache, paste0("staging-", Sys.getpid()))
assert(!dir.exists(staging_root), paste("Staging directory exists:", staging_root))
dir.create(staging_root, recursive = TRUE)
complete <- FALSE
on.exit(if (!complete && dir.exists(staging_root))
  unlink(staging_root, recursive = TRUE, force = TRUE), add = TRUE)

prepared <- lapply(requested, prepare_source)
names(prepared) <- requested
if ("cattle" %in% requested) build_cattle(prepared$cattle, staging_root)
if ("pig" %in% requested) build_pig(prepared$pig, staging_root)

# Transactional directory installation; roll back every requested dataset if
# any rename fails.
backups <- setNames(paste0(targets[requested], ".backup-", Sys.getpid()), requested)
installed <- setNames(rep(FALSE, length(requested)), requested)
backed_up <- setNames(rep(FALSE, length(requested)), requested)
rollback <- TRUE
on.exit(if (rollback) {
  for (kind in requested) if (installed[[kind]] && dir.exists(targets[[kind]]))
    unlink(targets[[kind]], recursive = TRUE, force = TRUE)
  for (kind in requested) if (backed_up[[kind]] && dir.exists(backups[[kind]]))
    file.rename(backups[[kind]], targets[[kind]])
}, add = TRUE)
for (kind in requested) if (dir.exists(targets[[kind]])) {
  assert(force, paste("Refusing to replace", targets[[kind]]))
  stage_name <- if (kind == "cattle") "cattle_milk_lipids" else "pig_blood_lipids"
  existing_readme <- file.path(targets[[kind]], "README.md")
  if (file.exists(existing_readme)) {
    assert(file.copy(existing_readme,
                     file.path(staging_root, stage_name, "README.md"),
                     overwrite = FALSE),
           paste("Could not preserve", existing_readme))
  }
  assert(file.rename(targets[[kind]], backups[[kind]]),
         paste("Could not back up", targets[[kind]]))
  backed_up[[kind]] <- TRUE
}
for (kind in requested) {
  stage <- file.path(staging_root, if (kind == "cattle")
    "cattle_milk_lipids" else "pig_blood_lipids")
  assert(file.rename(stage, targets[[kind]]), paste("Could not install", targets[[kind]]))
  installed[[kind]] <- TRUE
}
for (kind in requested) if (backed_up[[kind]])
  unlink(backups[[kind]], recursive = TRUE, force = TRUE)
rollback <- FALSE
complete <- TRUE
unlink(staging_root, recursive = TRUE, force = TRUE)

cat("status=PASS\n")
cat("species=", paste(requested, collapse = ","), "\n", sep = "")
cat("plink2=", paste(run_plink2("--version", "version"), collapse = " "), "\n", sep = "")
for (kind in requested) {
  x <- prepared[[kind]]
  cat(kind, "_version_url=", x$url, "\n", sep = "")
  cat(kind, "_outer_bytes=", x$outer_bytes, "\n", sep = "")
  cat(kind, "_outer_sha256=", x$outer_sha256, "\n", sep = "")
  cat(kind, "_inner_bytes=", x$inner_bytes, "\n", sep = "")
  cat(kind, "_inner_", x$inner_digest_type, "=", x$inner_digest, "\n", sep = "")
}

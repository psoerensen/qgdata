# Shared, focused checkpoint and PLINK assembly helpers for the human panels.

qg_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

qg_sha256 <- function(path) {
  qg_assert(length(path) == 1L && file.exists(path),
            paste("Cannot hash missing file:", path))
  unname(digest::digest(path, algo = "sha256", file = TRUE))
}

qg_hash_values <- function(x) {
  digest::digest(enc2utf8(as.character(x)), algo = "sha256", serialize = TRUE)
}

qg_expected_bed_bytes <- function(n_individuals, n_markers) {
  3 + as.double(n_markers) * ceiling(as.double(n_individuals) / 4)
}

qg_valid_hard_calls <- function(x) {
  is.matrix(x) && !anyNA(x) && all(is.finite(x)) && all(x %in% 0:2)
}

qg_checkpoint_paths <- function(path) {
  c(simulated = file.path(path, "simulated.rds"),
    reference = file.path(path, "reference.rds"),
    metadata = file.path(path, "metadata.rds"))
}

qg_metadata_mismatches <- function(actual, expected) {
  missing <- setdiff(names(expected), names(actual))
  different <- setdiff(names(expected), missing)
  different <- different[!vapply(different, function(nm) {
    identical(actual[[nm]], expected[[nm]])
  }, logical(1))]
  missing_labels <- if (length(missing)) paste0(missing, " (missing)") else character()
  c(missing_labels, different)
}

qg_validate_checkpoint <- function(path, expected, read_matrices = TRUE) {
  paths <- qg_checkpoint_paths(path)
  qg_assert(dir.exists(path), paste("Checkpoint directory is missing:", path))
  qg_assert(all(file.exists(paths)), paste("Checkpoint is incomplete:", path))
  metadata <- readRDS(paths[["metadata"]])
  mismatch <- qg_metadata_mismatches(metadata, expected)
  qg_assert(!length(mismatch), paste0(
    "Checkpoint metadata mismatch at ", path, ": ",
    paste(mismatch, collapse = ", "),
    ". Use --force to regenerate it; --resume never reuses incompatible data."
  ))
  qg_assert(identical(metadata$simulated_sha256, qg_sha256(paths[["simulated"]])),
            paste("Checkpoint simulated-matrix checksum differs:", path))
  qg_assert(identical(metadata$reference_sha256, qg_sha256(paths[["reference"]])),
            paste("Checkpoint reference-matrix checksum differs:", path))
  if (!read_matrices) return(invisible(metadata))

  G <- readRDS(paths[["simulated"]])
  R <- readRDS(paths[["reference"]])
  expected_dim <- c(expected$n_individuals, expected$n_markers)
  qg_assert(identical(dim(G), expected_dim),
            paste("Checkpoint simulated dimensions differ:", path))
  qg_assert(identical(rownames(G), expected$individual_ids),
            paste("Checkpoint simulated individual IDs/order differ:", path))
  qg_assert(identical(colnames(G), expected$marker_ids),
            paste("Checkpoint simulated marker IDs/order differ:", path))
  qg_assert(qg_valid_hard_calls(G),
            paste("Checkpoint simulated data are not complete diploid 0/1/2 calls:", path))
  qg_assert(identical(ncol(R), expected$n_markers),
            paste("Checkpoint reference marker count differs:", path))
  qg_assert(identical(colnames(R), expected$marker_ids),
            paste("Checkpoint reference marker IDs/order differ:", path))
  qg_assert(identical(rownames(R), expected$reference_ids),
            paste("Checkpoint reference individual IDs/order differ:", path))
  qg_assert(qg_valid_hard_calls(R),
            paste("Checkpoint reference data are not complete diploid 0/1/2 calls:", path))
  invisible(list(metadata = metadata, simulated = G, reference = R))
}

qg_write_checkpoint_atomic <- function(path, G, R, metadata, replace = FALSE) {
  parent <- dirname(path)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(parent, paste0(basename(path), ".partial-", Sys.getpid()))
  qg_assert(!dir.exists(tmp), paste("Temporary checkpoint already exists:", tmp))
  dir.create(tmp, recursive = FALSE)
  completed <- FALSE
  on.exit(if (!completed && dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE),
          add = TRUE)
  paths <- qg_checkpoint_paths(tmp)
  saveRDS(G, paths[["simulated"]], compress = "gzip")
  saveRDS(R, paths[["reference"]], compress = "gzip")
  metadata$simulated_sha256 <- qg_sha256(paths[["simulated"]])
  metadata$reference_sha256 <- qg_sha256(paths[["reference"]])
  metadata$created_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  saveRDS(metadata, paths[["metadata"]], compress = "gzip")
  qg_validate_checkpoint(tmp, metadata[names(metadata) != "created_utc"])
  if (dir.exists(path)) {
    qg_assert(replace, paste("Checkpoint already exists:", path))
    unlink(path, recursive = TRUE, force = TRUE)
  }
  qg_assert(file.rename(tmp, path), paste("Atomic checkpoint rename failed:", path))
  completed <- TRUE
  invisible(metadata)
}

qg_plink_fam <- function(ids) {
  data.frame(fam = ids, id = ids, pat = "0", mat = "0", sex = 0L,
             pheno = -9, stringsAsFactors = FALSE)
}

qg_validate_plink_files <- function(prefix, selected, individual_ids) {
  files <- paste0(prefix, c(".bed", ".bim", ".fam"))
  qg_assert(all(file.exists(files)), paste("Incomplete staged PLINK fileset:", prefix))
  expected_size <- qg_expected_bed_bytes(length(individual_ids), nrow(selected))
  actual_size <- unname(file.info(paste0(prefix, ".bed"))$size)
  qg_assert(identical(as.double(actual_size), expected_size), paste0(
    "BED size differs for ", prefix, ": expected ", expected_size,
    " bytes, found ", actual_size
  ))
  bim <- genio::read_bim(prefix, verbose = FALSE)
  fam <- genio::read_fam(prefix, verbose = FALSE)
  qg_assert(nrow(bim) == nrow(selected), "BIM marker count differs")
  qg_assert(nrow(fam) == length(individual_ids), "FAM individual count differs")
  qg_assert(identical(as.character(bim$id), as.character(selected$ID)),
            "BIM marker IDs/order differ")
  qg_assert(identical(as.integer(bim$chr), as.integer(selected$CHROM)),
            "BIM chromosomes differ")
  qg_assert(identical(as.double(bim$pos), as.double(selected$POS)),
            "BIM physical positions differ")
  qg_assert(identical(as.character(fam$id), individual_ids),
            "FAM individual IDs/order differ")
  qg_assert(identical(as.character(fam$fam), individual_ids),
            "FAM family IDs/order differ")
  invisible(list(bed_bytes = expected_size, bim = bim, fam = fam))
}

qg_assemble_checkpoints <- function(checkpoint_paths, expected_metadata, selected,
                                    staging_parent, ancestry, use_qgg = TRUE) {
  qg_assert(length(checkpoint_paths) == length(expected_metadata),
            "Checkpoint and expected-metadata counts differ")
  qg_assert(identical(as.integer(sort(unique(selected$block))),
                      as.integer(seq_along(checkpoint_paths))),
            "Selected variants do not contain each expected block exactly once")
  qg_assert(!anyDuplicated(selected$ID), "Selected marker IDs are duplicated")
  individual_ids <- expected_metadata[[1]]$individual_ids
  qg_assert(all(vapply(expected_metadata, function(x)
    identical(x$individual_ids, individual_ids), logical(1))),
    "Individual IDs/order differ across blocks")

  dir.create(staging_parent, recursive = TRUE, showWarnings = FALSE)
  stage <- file.path(staging_parent, paste0(tolower(ancestry), ".partial-", Sys.getpid()))
  qg_assert(!dir.exists(stage), paste("Staging directory exists:", stage))
  dir.create(stage)
  complete <- FALSE
  on.exit(if (!complete && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE),
          add = TRUE)
  prefix <- file.path(stage, paste0("human_1000g_", tolower(ancestry)))
  fam <- qg_plink_fam(individual_ids)

  for (block in seq_along(checkpoint_paths)) {
    validated <- qg_validate_checkpoint(checkpoint_paths[[block]],
                                        expected_metadata[[block]])
    G <- validated$simulated
    meta <- selected[selected$block == block, , drop = FALSE]
    bim <- data.frame(chr = meta$CHROM, id = meta$ID, posg = meta$CM,
                      pos = meta$POS, ref = meta$REF, alt = meta$ALT,
                      stringsAsFactors = FALSE)
    genio::write_plink(prefix, t(G), bim = bim, fam = fam, verbose = FALSE,
                       append = block > 1L)
    rm(G, validated)
    gc(verbose = FALSE)
  }
  qg_validate_plink_files(prefix, selected, individual_ids)

  if (use_qgg && requireNamespace("qgg", quietly = TRUE)) {
    qgg_dir <- file.path(stage, "qgg")
    dir.create(qgg_dir)
    old <- setwd(qgg_dir)
    on.exit(setwd(old), add = TRUE)
    glist <- qgg::gprep(study = paste0("assemble_", tolower(ancestry)),
                        bedfiles = paste0(prefix, ".bed"),
                        bimfiles = paste0(prefix, ".bim"),
                        famfiles = paste0(prefix, ".fam"), overwrite = TRUE)
    setwd(old)
    for (block in seq_along(checkpoint_paths)) {
      cls <- which(selected$block == block)
      back <- qgg::getG(Glist = glist, chr = 1, cls = cls,
                        impute = FALSE, scale = FALSE)
      source <- readRDS(qg_checkpoint_paths(checkpoint_paths[[block]])[["simulated"]])
      qg_assert(identical(dim(back), dim(source)), "qgg dimensions differ")
      qg_assert(identical(rownames(back), rownames(source)), "qgg individual order differs")
      qg_assert(identical(colnames(back), colnames(source)), "qgg marker order differs")
      qg_assert(identical(is.na(back), is.na(source)), "qgg missingness differs")
      qg_assert(isTRUE(all.equal(unname(back), unname(source), tolerance = 0,
                                check.attributes = FALSE)),
                "qgg genotype values differ")
    }
  }
  complete <- TRUE
  list(directory = stage, prefix = prefix,
       bed_bytes = qg_expected_bed_bytes(length(individual_ids), nrow(selected)))
}

qg_install_plink_transaction <- function(staged_prefix, destination_prefix,
                                         force = FALSE, fail_before_install = FALSE) {
  extensions <- c(".bed", ".bim", ".fam")
  source <- paste0(staged_prefix, extensions)
  destination <- paste0(destination_prefix, extensions)
  qg_assert(all(file.exists(source)), "Staged PLINK files are incomplete")
  if (fail_before_install) stop("Intentional test failure before canonical installation",
                                call. = FALSE)
  qg_assert(force || !any(file.exists(destination)), paste0(
    "Canonical PLINK output already exists at ", destination_prefix,
    "; use --force only when intentional"
  ))
  dir.create(dirname(destination_prefix), recursive = TRUE, showWarnings = FALSE)
  token <- paste0(".installing-", Sys.getpid())
  temporary <- paste0(destination_prefix, token, extensions)
  backup <- paste0(destination_prefix, ".backup-", Sys.getpid(), extensions)
  qg_assert(!any(file.exists(c(temporary, backup))), "Stale install transaction files exist")
  cleanup <- TRUE
  installed <- rep(FALSE, length(extensions))
  backed_up <- rep(FALSE, length(extensions))
  on.exit({
    if (cleanup) {
      unlink(temporary, force = TRUE)
      unlink(destination[installed], force = TRUE)
      for (i in which(backed_up)) if (file.exists(backup[[i]]))
        file.rename(backup[[i]], destination[[i]])
    }
  }, add = TRUE)
  for (i in seq_along(source)) {
    qg_assert(file.copy(source[[i]], temporary[[i]], overwrite = FALSE),
              paste("Could not stage canonical install:", temporary[[i]]))
    qg_assert(identical(qg_sha256(source[[i]]), qg_sha256(temporary[[i]])),
              paste("Install copy checksum differs:", temporary[[i]]))
  }
  for (i in seq_along(destination)) if (file.exists(destination[[i]])) {
    qg_assert(file.rename(destination[[i]], backup[[i]]),
              paste("Could not back up existing canonical file:", destination[[i]]))
    backed_up[[i]] <- TRUE
  }
  for (i in seq_along(destination)) {
    qg_assert(file.rename(temporary[[i]], destination[[i]]),
              paste("Could not install canonical file:", destination[[i]]))
    installed[[i]] <- TRUE
  }
  qg_assert(all(vapply(seq_along(source), function(i)
    identical(qg_sha256(source[[i]]), qg_sha256(destination[[i]])), logical(1))),
    "Installed PLINK checksum differs from staging")
  unlink(backup, force = TRUE)
  cleanup <- FALSE
  invisible(destination)
}

qg_atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".partial-", Sys.getpid())
  backup <- paste0(path, ".backup-", Sys.getpid())
  restore <- FALSE
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
    if (restore && file.exists(backup)) file.rename(backup, path)
  }, add = TRUE)
  data.table::fwrite(x, tmp)
  if (file.exists(path)) {
    qg_assert(file.rename(path, backup), paste("Could not back up table:", path))
    restore <- TRUE
  }
  qg_assert(file.rename(tmp, path), paste("Atomic table rename failed:", path))
  if (file.exists(backup)) unlink(backup)
  restore <- FALSE
  invisible(path)
}

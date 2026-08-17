#!/usr/bin/env Rscript

# Locked ancestry-specific R haplotype-mosaic scaling and production build.
# Simulation is blockwise and marker-major, but implements the same donor-copy
# transition as simulate_haplotype_mosaic(): 1-exp(-rate * delta_morgan).

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}
mode <- match.arg(tolower(arg_value("mode", "scaling")), c("scaling", "full"))
ancestry_arg <- match.arg(toupper(arg_value("ancestry", "both")),
                          c("EUR", "AFR", "BOTH"))
ancestries <- if (ancestry_arg == "BOTH") c("EUR", "AFR") else ancestry_arg
resume <- "--resume" %in% args
force <- "--force" %in% args
if (resume && force) stop("--resume and --force are mutually exclusive")
if (mode == "full" && !identical(ancestries, c("EUR", "AFR")))
  stop("Production installation requires --ancestry=both so both panels pass before installation")

required <- c("data.table", "digest", "genio", "qgg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "qgdata.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find qgdata repository root")
    path <- parent
  }
}
root <- find_root()
source(file.path(root, "scripts", "human_1000g_common.R"), local = TRUE)
source(file.path(root, "scripts", "human_panel_validation_common.R"), local = TRUE)
source(file.path(root, "scripts", "simulate_haplotype_mosaic.R"), local = TRUE)

cache_1000g <- file.path(root, "cache", "1000g_phase3")
build_root <- file.path(root, "cache", "human_r_mosaic",
                        if (mode == "full") "production" else "scaling")
validation_dir <- file.path(root, "validation")
evidence_prefix <- if (mode == "full") "human_validation" else "human_mosaic_scaling"
evidence_path <- function(kind) file.path(validation_dir,
  paste0(evidence_prefix, "_", kind, ".csv"))
plink2 <- file.path(root, "cache", "tools", "plink2_20260808", "plink2.exe")
map_path <- file.path(cache_1000g, "genetic_map_GRCh37_chr22.txt.gz")
source_prefix <- file.path(cache_1000g, "source", "chr22_phase3_biallelic_snps")
qg_assert(file.exists(plink2), paste("Missing cached PLINK2:", plink2))
qg_assert(file.exists(map_path), paste("Missing cached GRCh37 map:", map_path))
qg_assert(all(file.exists(paste0(source_prefix, c(".pgen", ".pvar", ".psam")))),
          paste("Missing cached phased PGEN source:", source_prefix))
dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, showWarnings = FALSE)

rates <- c(EUR = 20, AFR = 25)
base_seed <- 20260817L
build_seed <- function(mode, ancestry) base_seed +
  (if (mode == "scaling") 600000L else 700000L) +
  (if (ancestry == "EUR") 1L else 2L)
validation_seed <- if (mode == "scaling") 202608177L else 202608178L
n_individuals <- if (mode == "scaling") 1000L else 5000L
n_markers <- if (mode == "scaling") 10000L else 50000L
block_size <- if (mode == "scaling") 10000L else 2500L
n_blocks <- as.integer(n_markers / block_size)
qg_assert(n_blocks * block_size == n_markers, "Marker count is not block-divisible")

run_plink2 <- function(arguments, label) {
  output <- system2(plink2, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L)
    stop("PLINK2 failed during ", label, ":\n", paste(output, collapse = "\n"))
  invisible(output)
}
plink_version <- paste(run_plink2("--version", "version check"), collapse = " ")

safe_cache_remove <- function(path) {
  root_resolved <- normalizePath(build_root, winslash = "/", mustWork = TRUE)
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  qg_assert(startsWith(paste0(target, "/"), paste0(root_resolved, "/")),
            paste("Refusing to remove outside build cache:", target))
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  else if (file.exists(path)) unlink(path, force = TRUE)
}

map <- data.table::fread(map_path, data.table = FALSE)
map_bp <- as.numeric(map[[2L]]); map_cm <- as.numeric(map[[4L]])
if (mode == "scaling") {
  selected_path <- file.path(root, "cache", "human_r_mosaic",
    "relatedness_refinement", "diagnostic_10k", "selected_variants.csv")
  qg_assert(file.exists(selected_path), paste("Missing locked 10,000-marker design:",
                                             selected_path))
  selected <- data.table::fread(selected_path, data.table = FALSE)
  selected$block <- 1L
  design_seed <- 202608175L
  marker_design <- "locked shared-MAF natural-density 10,000-marker diagnostic design"
} else {
  selected_path <- file.path(cache_1000g, "build_full", "selected_variants.csv")
  qg_assert(file.exists(selected_path), paste("Missing locked 50,000-marker design:",
                                             selected_path))
  selected <- data.table::fread(selected_path, data.table = FALSE)
  selected$CM <- stats::approx(map_bp, map_cm, selected$POS, rule = 2)$y
  design_seed <- base_seed
  marker_design <- "documented deterministic broad-coverage 50,000-marker design"
}
qg_assert(nrow(selected) == n_markers && !anyDuplicated(selected$ID) &&
  all(selected$CHROM == 22L) && all(diff(selected$POS) >= 0) &&
  all(is.finite(selected$CM)) && all(diff(selected$CM) >= 0),
  "Locked marker design is invalid")
selected$block <- ceiling(seq_len(nrow(selected)) / block_size)
marker_checksum <- qg_hash_values(selected$ID)
source_checksums <- c(pgen = qg_sha256(paste0(source_prefix, ".pgen")),
  pvar = qg_sha256(paste0(source_prefix, ".pvar")),
  psam = qg_sha256(paste0(source_prefix, ".psam")),
  map = qg_sha256(map_path), selected = qg_sha256(selected_path))

read_phased_vcf_fast <- function(path, meta, expected_samples) {
  qg_assert(file.exists(path), paste("Missing phased reference:", path))
  connection <- gzfile(path, "rt")
  header <- readLines(connection, n = 500L, warn = FALSE)
  close(connection)
  header_line <- grep("^#CHROM\\t", header, value = TRUE)
  qg_assert(length(header_line) == 1L, paste("Invalid VCF header:", path))
  header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  samples <- header_fields[-seq_len(9L)]
  qg_assert(length(samples) == expected_samples && !anyDuplicated(samples),
            paste("Unexpected reference sample count:", path))
  x <- data.table::fread(path, skip = "#CHROM", data.table = FALSE,
                         showProgress = FALSE)
  names(x)[names(x) == "#CHROM"] <- "CHROM"
  qg_assert(nrow(x) == nrow(meta) &&
    identical(as.character(x$ID), as.character(meta$ID)) &&
    identical(as.integer(x$POS), as.integer(meta$POS)) &&
    identical(as.character(x$REF), as.character(meta$REF)) &&
    identical(as.character(x$ALT), as.character(meta$ALT)),
    paste("VCF marker identity/order/alleles differ:", path))
  qg_assert(all(x$FORMAT == "GT"), paste("Reference VCF FORMAT is not exactly GT:", path))
  gt <- as.matrix(x[, samples, drop = FALSE])
  valid <- nchar(gt) == 3L & substr(gt, 2L, 2L) == "|" &
    substr(gt, 1L, 1L) %in% c("0", "1") & substr(gt, 3L, 3L) %in% c("0", "1")
  qg_assert(all(valid), paste("Reference contains missing/nonbinary/unphased GT:", path))
  a1 <- matrix(as.integer(substr(gt, 1L, 1L)), nrow = nrow(gt), ncol = ncol(gt))
  a2 <- matrix(as.integer(substr(gt, 3L, 3L)), nrow = nrow(gt), ncol = ncol(gt))
  H <- matrix(0L, 2L * length(samples), nrow(meta),
    dimnames = list(as.vector(rbind(paste0(samples, "_H1"), paste0(samples, "_H2"))),
                    meta$ID))
  H[seq.int(1L, nrow(H), 2L), ] <- t(a1)
  H[seq.int(2L, nrow(H), 2L), ] <- t(a2)
  R <- H[seq.int(1L, nrow(H), 2L), , drop = FALSE] +
    H[seq.int(2L, nrow(H), 2L), , drop = FALSE]
  rownames(R) <- samples; colnames(R) <- meta$ID; storage.mode(R) <- "integer"
  list(H = H, R = R, samples = samples, sha256 = qg_sha256(path))
}

prepare_reference_block <- function(ancestry, block, meta) {
  reference_dir <- file.path(build_root, "reference", tolower(ancestry))
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(reference_dir, sprintf("block_%02d", block))
  vcf <- paste0(prefix, ".vcf.gz")
  if (mode == "scaling") {
    existing <- file.path(root, "cache", "human_r_mosaic", "relatedness_refinement",
      "diagnostic_10k", paste0(tolower(ancestry), "_reference.vcf.gz"))
    qg_assert(file.exists(existing), paste("Missing locked scaling reference:", existing))
    vcf <- existing
  } else if (!file.exists(vcf)) {
    ids <- paste0(prefix, ".ids")
    data.table::fwrite(data.frame(ID = meta$ID), ids, col.names = FALSE)
    run_plink2(c("--pfile", source_prefix, "--keep",
      file.path(cache_1000g, "build_full", paste0(tolower(ancestry), ".keep")),
      "--extract", ids, "--export", "vcf", "bgz", "--out", prefix),
      paste(ancestry, "reference block", block))
  }
  read_phased_vcf_fast(vcf, meta, if (ancestry == "EUR") 503L else 504L)
}

checkpoint_paths <- function(directory) c(genotypes = file.path(directory, "genotypes.rds"),
  reference = file.path(directory, "reference.rds"), state = file.path(directory, "state.rds"),
  metadata = file.path(directory, "metadata.rds"),
  donors = file.path(directory, "donors.rds"))

write_block_atomic <- function(directory, G, R, state, metadata, donors = NULL) {
  temporary <- paste0(directory, ".tmp-", Sys.getpid())
  qg_assert(!dir.exists(temporary), paste("Temporary checkpoint exists:", temporary))
  dir.create(temporary, recursive = TRUE)
  complete <- FALSE
  on.exit(if (!complete && dir.exists(temporary))
    unlink(temporary, recursive = TRUE, force = TRUE), add = TRUE)
  saveRDS(G, file.path(temporary, "genotypes.rds"), compress = FALSE)
  saveRDS(R, file.path(temporary, "reference.rds"), compress = FALSE)
  saveRDS(state, file.path(temporary, "state.rds"), compress = FALSE)
  saveRDS(metadata, file.path(temporary, "metadata.rds"))
  if (!is.null(donors)) saveRDS(donors, file.path(temporary, "donors.rds"), compress = FALSE)
  check <- readRDS(file.path(temporary, "genotypes.rds"))
  qg_assert(identical(dim(check), c(metadata$n_individuals, metadata$n_markers)) &&
    identical(rownames(check), metadata$individual_ids) &&
    identical(colnames(check), metadata$marker_ids) && qg_valid_hard_calls(check),
    "Temporary genotype checkpoint failed validation")
  qg_assert(file.rename(temporary, directory), paste("Could not finalize checkpoint:", directory))
  complete <- TRUE
}

validate_checkpoint <- function(directory, expected) {
  paths <- checkpoint_paths(directory)
  needed <- paths[names(paths) != "donors"]
  qg_assert(all(file.exists(needed)), paste("Incomplete block checkpoint:", directory))
  metadata <- readRDS(paths[["metadata"]])
  qg_assert(identical(metadata, expected), paste("Incompatible block checkpoint:", directory))
  G <- readRDS(paths[["genotypes"]]); R <- readRDS(paths[["reference"]])
  qg_assert(identical(dim(G), c(expected$n_individuals, expected$n_markers)) &&
    identical(rownames(G), expected$individual_ids) &&
    identical(colnames(G), expected$marker_ids) && qg_valid_hard_calls(G),
    paste("Invalid simulated checkpoint:", directory))
  qg_assert(identical(dim(R), c(expected$n_reference, expected$n_markers)) &&
    identical(rownames(R), expected$reference_ids) &&
    identical(colnames(R), expected$marker_ids) && qg_valid_hard_calls(R),
    paste("Invalid reference checkpoint:", directory))
  list(G = G, R = R, state = readRDS(paths[["state"]]), metadata = metadata,
       donors = if (file.exists(paths[["donors"]])) readRDS(paths[["donors"]]) else NULL)
}

build_ancestry <- function(ancestry) {
  ancestry_dir <- file.path(build_root, tolower(ancestry))
  if (dir.exists(ancestry_dir) && force) safe_cache_remove(ancestry_dir)
  if (dir.exists(ancestry_dir) && !resume)
    stop("Build cache already exists for ", ancestry, "; use --resume or --force")
  dir.create(ancestry_dir, recursive = TRUE, showWarnings = FALSE)
  individual_ids <- sprintf("MOSAIC_%s_%05d", ancestry, seq_len(n_individuals))
  seed <- build_seed(mode, ancestry)
  state <- NULL
  block_seconds <- numeric(n_blocks)
  reference_hashes <- character(n_blocks)
  donor_blocks <- if (mode == "scaling") vector("list", n_blocks) else NULL
  for (block in seq_len(n_blocks)) {
    meta <- selected[selected$block == block, , drop = FALSE]
    reference <- prepare_reference_block(ancestry, block, meta)
    reference_hashes[[block]] <- reference$sha256
    expected <- list(schema_version = 1L, mode = mode, ancestry = ancestry,
      copy_rate = unname(rates[[ancestry]]), seed = seed, block = as.integer(block),
      n_individuals = n_individuals, n_markers = nrow(meta),
      n_reference = nrow(reference$R), individual_ids = individual_ids,
      reference_ids = rownames(reference$R), marker_ids = as.character(meta$ID),
      marker_checksum = qg_hash_values(meta$ID), full_marker_checksum = marker_checksum,
      reference_sha256 = reference$sha256, source_pgen_sha256 = source_checksums[["pgen"]],
      genetic_map_sha256 = source_checksums[["map"]], genome_build = "GRCh37",
      chromosome = 22L, block_start = min(which(selected$block == block)),
      block_end = max(which(selected$block == block)),
      block_start_bp = min(meta$POS), block_end_bp = max(meta$POS))
    checkpoint <- file.path(ancestry_dir, "blocks", sprintf("block_%02d", block))
    if (dir.exists(checkpoint)) {
      qg_assert(resume, paste("Checkpoint exists without --resume:", checkpoint))
      reused <- validate_checkpoint(checkpoint, expected)
      state <- reused$state
      if (mode == "scaling") donor_blocks[[block]] <- reused$donors
      message("REUSE ", ancestry, " block ", block)
      next
    }
    if (is.null(state)) {
      state <- initialize_haplotype_mosaic_state(2L * n_individuals,
        nrow(reference$H), seed, meta$CM[[1L]] / 100, meta$POS[[1L]])
      state$reference_haplotype_ids <- rownames(reference$H)
    } else {
      qg_assert(state$n_donors == nrow(reference$H), "Reference donor count changed")
      qg_assert(identical(state$reference_haplotype_ids, rownames(reference$H)),
                "Reference haplotype identity/order changed across blocks")
      assign(".Random.seed", state$rng_state, envir = .GlobalEnv)
    }
    start <- proc.time()[["elapsed"]]
    simulated <- simulate_haplotype_mosaic_block(reference$H, meta$CM / 100,
      meta$POS, unname(rates[[ancestry]]), state, return_donors = mode == "scaling")
    block_seconds[[block]] <- proc.time()[["elapsed"]] - start
    G <- simulated$genotypes
    rownames(G) <- individual_ids; colnames(G) <- meta$ID
    state <- simulated$state; state$rng_state <- get(".Random.seed", envir = .GlobalEnv)
    donors <- simulated$donors
    write_block_atomic(checkpoint, G, reference$R, state, expected, donors)
    if (mode == "scaling") donor_blocks[[block]] <- donors
    message("GENERATE ", ancestry, " block ", block,
            sprintf(" (%.2fs)", block_seconds[[block]]))
    rm(G, simulated, reference); gc(verbose = FALSE)
  }
  list(directory = ancestry_dir, individual_ids = individual_ids, seed = seed,
       reference_hashes = reference_hashes, block_seconds = block_seconds,
       donor_blocks = donor_blocks,
       segment_summary = finalize_haplotype_mosaic_state(state))
}

assemble_ancestry <- function(ancestry, build) {
  assembly <- file.path(build$directory, "assembly")
  if (dir.exists(assembly)) {
    qg_assert(resume, paste("Assembly exists without --resume:", assembly))
    safe_cache_remove(assembly)
  }
  dir.create(assembly, recursive = TRUE)
  prefix <- file.path(assembly, paste0("human_1000g_", tolower(ancestry)))
  fam <- qg_plink_fam(build$individual_ids)
  duplicate_fingerprint <- rep("", n_individuals)
  full_maf <- numeric(n_markers)
  heterozygous <- 0
  validation_index <- if (mode == "scaling") seq_len(n_markers) else
    unique(floor(seq(1, n_markers, length.out = 10000L)))
  qg_assert(length(validation_index) == 10000L, "Validation marker selection is not 10,000")
  G_validation <- matrix(0L, n_individuals, length(validation_index),
    dimnames = list(build$individual_ids, selected$ID[validation_index]))
  R_validation <- NULL
  validation_cursor <- 1L
  for (block in seq_len(n_blocks)) {
    checkpoint <- file.path(build$directory, "blocks", sprintf("block_%02d", block))
    paths <- checkpoint_paths(checkpoint)
    G <- readRDS(paths[["genotypes"]]); R <- readRDS(paths[["reference"]])
    meta <- selected[selected$block == block, , drop = FALSE]
    block_global <- which(selected$block == block)
    full_maf[block_global] <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
    heterozygous <- heterozygous + sum(G == 1L)
    duplicate_fingerprint <- paste0(duplicate_fingerprint, qg_row_hashes(G))
    take_global <- intersect(validation_index, block_global)
    if (length(take_global)) {
      local <- match(take_global, block_global)
      destination <- match(take_global, validation_index)
      G_validation[, destination] <- G[, local, drop = FALSE]
      if (is.null(R_validation)) R_validation <- matrix(0L, nrow(R), length(validation_index),
        dimnames = list(rownames(R), selected$ID[validation_index]))
      R_validation[, destination] <- R[, local, drop = FALSE]
    }
    bim <- data.frame(chr = meta$CHROM, id = meta$ID, posg = meta$CM,
      pos = meta$POS, ref = meta$REF, alt = meta$ALT, stringsAsFactors = FALSE)
    genio::write_plink(prefix, t(G), bim = bim, fam = fam, verbose = FALSE,
                       append = block > 1L)
    rm(G, R); gc(verbose = FALSE)
  }
  qg_validate_plink_files(prefix, selected, build$individual_ids)
  expected_size <- qg_expected_bed_bytes(n_individuals, n_markers)
  qg_assert(file.info(paste0(prefix, ".bed"))$size == expected_size,
            "Final staged BED size differs from exact expectation")
  list(prefix = prefix, full_maf = full_maf,
    heterozygosity = heterozygous / (n_individuals * n_markers),
    duplicated_individuals = sum(duplicated(vapply(duplicate_fingerprint,
      digest::digest, character(1), algo = "sha256", serialize = FALSE))),
    G_validation = G_validation, R_validation = R_validation,
    validation_index = validation_index, expected_bed_bytes = expected_size)
}

king_summary <- function(prefix, n, label) {
  out <- file.path(dirname(prefix), paste0(basename(prefix), "_king"))
  run_plink2(c("--bfile", prefix, "--make-king", "triangle", "bin4", "--out", out),
             paste(label, "KING"))
  expected_pairs <- n * (n - 1) / 2
  qg_assert(file.info(paste0(out, ".king.bin"))$size == expected_pairs * 4,
            paste(label, "KING binary size differs"))
  values <- readBin(paste0(out, ".king.bin"), what = "numeric", n = expected_pairs,
                    size = 4L)
  qg_assert(length(values) == expected_pairs && all(is.finite(values)),
            paste(label, "KING coefficients are incomplete"))
  c(pairs = expected_pairs, median = stats::median(values),
    q950 = unname(stats::quantile(values, .95)),
    q975 = unname(stats::quantile(values, .975)),
    q995 = unname(stats::quantile(values, .995)), maximum = max(values),
    duplicate_MZ_pairs = sum(values >= 2^(-1.5)),
    first_degree_pairs = sum(values >= 2^(-2.5) & values < 2^(-1.5)),
    second_degree_pairs = sum(values >= 2^(-3.5) & values < 2^(-2.5)))
}

focused_roundtrip <- function(ancestry, assembled) {
  focus <- assembled$validation_index[unique(round(seq(1,
    length(assembled$validation_index), length.out = 100L)))]
  ids <- file.path(dirname(assembled$prefix), "focused.ids")
  data.table::fwrite(data.frame(ID = selected$ID[focus]), ids, col.names = FALSE)
  out <- file.path(dirname(assembled$prefix), "focused_roundtrip")
  run_plink2(c("--bfile", assembled$prefix, "--extract", ids, "--make-bed", "--out", out),
             paste(ancestry, "focused extraction"))
  back <- t(genio::read_plink(out, verbose = FALSE)$X)
  expected <- assembled$G_validation[, match(focus, assembled$validation_index), drop = FALSE]
  qg_assert(identical(dim(back), dim(expected)) &&
    identical(rownames(back), rownames(expected)) &&
    identical(colnames(back), colnames(expected)) &&
    isTRUE(all.equal(unname(back), unname(expected), tolerance = 0,
                     check.attributes = FALSE)), "Focused genio round-trip differs")
  qgg_dir <- file.path(dirname(assembled$prefix), "qgg")
  dir.create(qgg_dir, showWarnings = FALSE)
  old <- setwd(qgg_dir); on.exit(setwd(old), add = TRUE)
  glist <- qgg::gprep(study = paste0("production_", tolower(ancestry)),
    bedfiles = paste0(assembled$prefix, ".bed"),
    bimfiles = paste0(assembled$prefix, ".bim"),
    famfiles = paste0(assembled$prefix, ".fam"), overwrite = TRUE)
  setwd(old)
  qgg_back <- qgg::getG(Glist = glist, chr = 1, cls = focus,
                        impute = FALSE, scale = FALSE)
  qg_assert(identical(dim(qgg_back), dim(expected)) &&
    identical(rownames(qgg_back), rownames(expected)) &&
    identical(colnames(qgg_back), colnames(expected)) &&
    isTRUE(all.equal(unname(qgg_back), unname(expected), tolerance = 0,
                     check.attributes = FALSE)), "Focused qgg round-trip differs")
  TRUE
}

scientific_validation <- function(ancestry, assembled, full_king) {
  G <- assembled$G_validation; R <- assembled$R_validation
  positions <- selected$POS[assembled$validation_index]
  full_reference_maf <- pmin(colMeans(R) / 2, 1 - colMeans(R) / 2)
  pca <- qg_prepare_pca_reference(R)
  ref_indices <- qg_matched_reference_indices(nrow(R), 20L, 200L,
    validation_seed + if (ancestry == "EUR") 100001L else 200001L)
  set.seed(validation_seed + if (ancestry == "EUR") 300001L else 400001L)
  synthetic_indices <- replicate(5L, sort(sample.int(nrow(G), 200L)), simplify = FALSE)
  reference_metrics <- do.call(rbind, lapply(ref_indices, function(index)
    qg_dataset_metrics_scalable(R[index, , drop = FALSE], positions,
                                full_reference_maf, pca)))
  synthetic_metrics <- do.call(rbind, lapply(synthetic_indices, function(index)
    qg_dataset_metrics_scalable(G[index, , drop = FALSE], positions,
                                full_reference_maf, pca)))
  ref_related <- do.call(rbind, lapply(ref_indices, function(index)
    qg_kinship_distribution(R[index, , drop = FALSE])))
  syn_related <- do.call(rbind, lapply(synthetic_indices, function(index)
    qg_kinship_distribution(G[index, , drop = FALSE])))
  reference_maf_all <- if (ancestry == "EUR") selected$MAF_EUR else selected$MAF_AFR
  qg_assert(length(reference_maf_all) == n_markers && all(is.finite(reference_maf_all)),
            "Full ancestry reference MAF vector is unavailable")
  maf_correlation <- stats::cor(assembled$full_maf, reference_maf_all)
  mean_maf_difference <- abs(mean(assembled$full_maf) - mean(reference_maf_all))
  reference_heterozygosity <- mean(R == 1L)
  validation_heterozygosity <- mean(G == 1L)
  heterozygosity_difference <- abs(validation_heterozygosity - reference_heterozygosity)
  ratio <- function(metric) synthetic_metrics[, metric] / mean(reference_metrics[, metric])
  ld_names <- paste0("mean_r2_", qg_validation_distance_names[1:3])
  ld_ratios <- sapply(ld_names, ratio)
  eigen_effective <- ratio("effective_rank")
  eigen_largest <- ratio("largest_eigenvalue")
  ld_score <- ratio("median_ld_score")
  reference_first_rate_q975 <- unname(stats::quantile(
    ref_related[, "first_degree_pairs"] / ref_related[, "pairs"], .975))
  reference_mz_rate_q975 <- unname(stats::quantile(
    ref_related[, "duplicate_MZ_pairs"] / ref_related[, "pairs"], .975))
  q995_limit <- unname(stats::quantile(ref_related[, "q995"], .975))
  full_first_rate <- full_king[["first_degree_pairs"]] / full_king[["pairs"]]
  full_mz_rate <- full_king[["duplicate_MZ_pairs"]] / full_king[["pairs"]]
  gates <- data.frame(gate = c("dimensions_hard_calls", "marker_integrity",
    "maf_correlation", "mean_maf_difference", "heterozygosity_difference",
    "ld_below_100kb", "ld_score", "effective_rank", "largest_eigenvalue",
    "pca_displacement", "duplicate_individuals", "empirical_relatedness",
    "plink_roundtrip"), pass = c(
      identical(dim(G), c(n_individuals, 10000L)) && qg_valid_hard_calls(G),
      identical(colnames(G), selected$ID[assembled$validation_index]),
      maf_correlation >= .98, mean_maf_difference <= .01,
      heterozygosity_difference <= .01,
      all(is.finite(ld_ratios) & ld_ratios >= .8 & ld_ratios <= 1.2),
      all(is.finite(ld_score) & ld_score >= .8 & ld_score <= 1.2),
      all(is.finite(eigen_effective) & eigen_effective >= .8 & eigen_effective <= 1.2),
      all(is.finite(eigen_largest) & eigen_largest >= .8 & eigen_largest <= 1.2),
      all(synthetic_metrics[, "pca_centroid_rms_z"] <= .25),
      assembled$duplicated_individuals == 0,
      full_first_rate <= reference_first_rate_q975 &&
        full_mz_rate <= reference_mz_rate_q975,
      TRUE), stringsAsFactors = FALSE)
  summary <- data.frame(ancestry = ancestry, mode = mode,
    individuals = n_individuals, markers = n_markers,
    validation_individuals = 200L, validation_markers = 10000L,
    copy_rate = unname(rates[[ancestry]]), build_seed = build_seed(mode, ancestry),
    validation_seed = validation_seed, maf_correlation = maf_correlation,
    mean_maf_difference = mean_maf_difference,
    heterozygosity = assembled$heterozygosity,
    validation_heterozygosity = validation_heterozygosity,
    reference_heterozygosity = reference_heterozygosity,
    heterozygosity_difference = heterozygosity_difference,
    mean_r2_ratio_0_10kb_median = stats::median(ld_ratios[, 1L]),
    mean_r2_ratio_0_10kb_min = min(ld_ratios[, 1L]),
    mean_r2_ratio_0_10kb_max = max(ld_ratios[, 1L]),
    mean_r2_ratio_10_50kb_median = stats::median(ld_ratios[, 2L]),
    mean_r2_ratio_10_50kb_min = min(ld_ratios[, 2L]),
    mean_r2_ratio_10_50kb_max = max(ld_ratios[, 2L]),
    mean_r2_ratio_50_100kb_median = stats::median(ld_ratios[, 3L]),
    mean_r2_ratio_50_100kb_min = min(ld_ratios[, 3L]),
    mean_r2_ratio_50_100kb_max = max(ld_ratios[, 3L]),
    median_ld_score_ratio = stats::median(ld_score),
    effective_rank_ratio = stats::median(eigen_effective),
    largest_eigenvalue_ratio = stats::median(eigen_largest),
    pca_rms_z_max = max(synthetic_metrics[, "pca_centroid_rms_z"]),
    duplicated_individuals = assembled$duplicated_individuals,
    full_king_q995 = full_king[["q995"]],
    synthetic_subset_q995_median = stats::median(syn_related[, "q995"]),
    matched_reference_q995_q975 = q995_limit,
    q995_descriptive_warning = stats::median(syn_related[, "q995"]) > q995_limit,
    full_MZ_pairs = full_king[["duplicate_MZ_pairs"]],
    full_first_degree_pairs = full_king[["first_degree_pairs"]],
    full_second_degree_pairs = full_king[["second_degree_pairs"]],
    full_MZ_rate = full_mz_rate, matched_reference_MZ_rate_q975 = reference_mz_rate_q975,
    full_first_degree_rate = full_first_rate,
    matched_reference_first_degree_rate_q975 = reference_first_rate_q975,
    bed_bytes = file.info(paste0(assembled$prefix, ".bed"))$size,
    expected_bed_bytes = assembled$expected_bed_bytes,
    all_gates_pass = all(gates$pass), stringsAsFactors = FALSE)
  list(summary = summary, gates = transform(gates, ancestry = ancestry, mode = mode),
       reference_metrics = reference_metrics, synthetic_metrics = synthetic_metrics,
       reference_relatedness = ref_related, synthetic_relatedness = syn_related)
}

sampled_donor_summary <- function(build, assembled, ancestry) {
  if (mode != "scaling") return(NULL)
  donor <- do.call(cbind, build$donor_blocks)
  qg_assert(identical(dim(donor), c(2L * n_individuals, n_markers)),
            "Scaling donor matrix dimensions differ")
  data.frame(ancestry = ancestry,
    sampled_individual_donor_overlap_summary(donor, selected$CM / 100,
      selected$POS, 2000L,
      validation_seed + if (ancestry == "EUR") 500001L else 600001L),
    stringsAsFactors = FALSE)
}

install_both_atomically <- function(results) {
  extensions <- c(".bed", ".bim", ".fam")
  destinations <- unlist(lapply(c("EUR", "AFR"), function(ancestry) {
    base <- file.path(root, paste0("human_1000g_", tolower(ancestry)),
      paste0("human_1000g_", tolower(ancestry)))
    paste0(base, extensions)
  }))
  qg_assert(!any(file.exists(destinations)),
            "Canonical human PLINK files already exist; refusing to overwrite")
  sources <- unlist(lapply(c("EUR", "AFR"), function(ancestry)
    paste0(results[[ancestry]]$assembled$prefix, extensions)))
  qg_assert(all(file.exists(sources)), "Staged production files are incomplete")
  temporary <- paste0(destinations, ".installing-", Sys.getpid())
  installed <- rep(FALSE, length(destinations)); moved <- rep(FALSE, length(sources))
  complete <- FALSE
  on.exit(if (!complete) {
    unlink(destinations[installed], force = TRUE)
    for (i in which(moved)) if (file.exists(temporary[[i]]))
      file.rename(temporary[[i]], sources[[i]])
  }, add = TRUE)
  for (i in seq_along(sources)) {
    dir.create(dirname(destinations[[i]]), recursive = TRUE, showWarnings = FALSE)
    qg_assert(file.rename(sources[[i]], temporary[[i]]), "Could not stage canonical install")
    moved[[i]] <- TRUE
  }
  for (i in seq_along(temporary)) {
    qg_assert(file.rename(temporary[[i]], destinations[[i]]), "Canonical install failed")
    installed[[i]] <- TRUE
  }
  complete <- TRUE
  destinations
}

start_time <- Sys.time()
results <- list(); summary_rows <- gates_rows <- donor_rows <- relatedness_rows <- list()
compact_relatedness <- function(ancestry, validation, king) {
  summarize <- function(x, source) data.frame(ancestry = ancestry, source = source,
    subsets = nrow(x), q950_mean = mean(x[, "q950"]),
    q975_mean = mean(x[, "q975"]), q995_mean = mean(x[, "q995"]),
    q995_median = stats::median(x[, "q995"]),
    q995_q025 = unname(stats::quantile(x[, "q995"], .025)),
    q995_q975 = unname(stats::quantile(x[, "q995"], .975)),
    maximum = max(x[, "maximum"]), MZ_pairs = sum(x[, "duplicate_MZ_pairs"]),
    first_degree_pairs = sum(x[, "first_degree_pairs"]),
    second_degree_pairs = sum(x[, "second_degree_pairs"]), stringsAsFactors = FALSE)
  rbind(summarize(validation$synthetic_relatedness, "synthetic_200_person_subsets"),
    summarize(validation$reference_relatedness, "matched_reference_200_person_subsets"),
    data.frame(ancestry = ancestry, source = paste0("full_", mode, "_panel_PLINK2"),
      subsets = 1L, q950_mean = NA, q975_mean = NA, q995_mean = king[["q995"]],
      q995_median = king[["q995"]], q995_q025 = NA, q995_q975 = NA,
      maximum = king[["maximum"]], MZ_pairs = king[["duplicate_MZ_pairs"]],
      first_degree_pairs = king[["first_degree_pairs"]],
      second_degree_pairs = king[["second_degree_pairs"]], stringsAsFactors = FALSE))
}
for (ancestry in ancestries) {
  message("BUILD ", mode, " ", ancestry, " rate=", rates[[ancestry]])
  build <- build_ancestry(ancestry)
  assembled <- assemble_ancestry(ancestry, build)
  focused_roundtrip(ancestry, assembled)
  king <- king_summary(assembled$prefix, n_individuals, paste(mode, ancestry))
  validation <- scientific_validation(ancestry, assembled, king)
  donor <- sampled_donor_summary(build, assembled, ancestry)
  summary_rows[[ancestry]] <- validation$summary
  gates_rows[[ancestry]] <- validation$gates
  if (!is.null(donor)) donor_rows[[ancestry]] <- donor
  relatedness_rows[[ancestry]] <- compact_relatedness(ancestry, validation, king)
  results[[ancestry]] <- list(build = build, assembled = assembled,
                              validation = validation, king = king)
  if (!all(validation$gates$pass)) {
    interim <- data.table::rbindlist(summary_rows, fill = TRUE)
    interim$elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    data.table::fwrite(interim, evidence_path("summary"))
    data.table::fwrite(data.table::rbindlist(gates_rows, fill = TRUE),
      evidence_path("acceptance"))
    if (length(donor_rows)) data.table::fwrite(data.table::rbindlist(donor_rows, fill = TRUE),
      evidence_path("donor_summary"))
    data.table::fwrite(data.table::rbindlist(relatedness_rows, fill = TRUE),
      evidence_path("relatedness"))
    stop(mode, " ", ancestry, " failed hard scientific gates; production not installed")
  }
}

summary_table <- data.table::rbindlist(summary_rows, fill = TRUE)
gates_table <- data.table::rbindlist(gates_rows, fill = TRUE)
donor_table <- data.table::rbindlist(donor_rows, fill = TRUE)
summary_table$elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
summary_table$peak_memory_bytes <- NA_real_
summary_table$plink2_version <- plink_version
summary_table$R_version <- R.version.string
summary_table$marker_design <- marker_design
summary_table$marker_design_seed <- design_seed
summary_table$marker_checksum <- marker_checksum
data.table::fwrite(summary_table, evidence_path("summary"))
data.table::fwrite(gates_table, evidence_path("acceptance"))
if (nrow(donor_table)) data.table::fwrite(donor_table,
  evidence_path("donor_summary"))
data.table::fwrite(data.table::rbindlist(relatedness_rows, fill = TRUE),
  evidence_path("relatedness"))

if (mode == "full") {
  qg_assert(all(summary_table$all_gates_pass), "Both production ancestries did not pass")
  installed <- install_both_atomically(results)
  # Validate the installed byte-for-byte files before declaring success.
  for (ancestry in c("EUR", "AFR")) {
    canonical <- file.path(root, paste0("human_1000g_", tolower(ancestry)),
      paste0("human_1000g_", tolower(ancestry)))
    qg_validate_plink_files(canonical, selected, results[[ancestry]]$build$individual_ids)
  }
  checksum_rows <- do.call(rbind, lapply(c("EUR", "AFR"), function(ancestry) {
    prefix <- file.path(root, paste0("human_1000g_", tolower(ancestry)),
      paste0("human_1000g_", tolower(ancestry)))
    files <- paste0(prefix, c(".bed", ".bim", ".fam"))
    data.frame(ancestry = ancestry, extension = c("bed", "bim", "fam"),
      path = substring(files, nchar(root) + 2L), bytes = file.info(files)$size,
      sha256 = vapply(files, qg_sha256, character(1)), stringsAsFactors = FALSE)
  }))
  data.table::fwrite(checksum_rows,
    file.path(validation_dir, "human_validation_checksums.csv"))
}

capture.output(sessionInfo(), file = file.path(validation_dir,
  paste0(evidence_prefix, "_sessionInfo.txt")))
cat("status=PASS mode=", mode, " ancestries=", paste(ancestries, collapse = ","), "\n",
    sep = "")
print(summary_table)

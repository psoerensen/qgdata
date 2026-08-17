#!/usr/bin/env Rscript

# Focused post-install validation for the canonical mosaic panels. Scientific
# metrics are produced during staged construction; this script verifies that
# the installed files still match that evidence and repeats small genio/qgg
# read-backs without loading the full genotype matrix.

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
required <- c("data.table", "digest", "genio", "qgg")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

source(file.path(root, "scripts", "human_1000g_common.R"), local = TRUE)

selected_path <- file.path(root, "cache", "1000g_phase3", "build_full",
                           "selected_variants.csv")
plink2 <- file.path(root, "cache", "tools", "plink2_20260808", "plink2.exe")
checksums_path <- file.path(root, "validation", "human_validation_checksums.csv")
summary_path <- file.path(root, "validation", "human_validation_summary.csv")
acceptance_path <- file.path(root, "validation", "human_validation_acceptance.csv")
qg_assert(file.exists(selected_path), paste("Missing locked marker design:", selected_path))
qg_assert(file.exists(plink2), paste("Missing cached PLINK2:", plink2))
qg_assert(file.exists(checksums_path) && file.exists(summary_path) &&
            file.exists(acceptance_path), "Missing compact production validation evidence")

selected <- data.table::fread(selected_path, data.table = FALSE)
qg_assert(nrow(selected) == 50000L && !anyDuplicated(selected$ID),
          "Locked marker design is not the expected unique 50,000 markers")
checksums <- data.table::fread(checksums_path, data.table = FALSE)
summary <- data.table::fread(summary_path, data.table = FALSE)
acceptance <- data.table::fread(acceptance_path, data.table = FALSE)
qg_assert(nrow(summary) == 2L && all(summary$individuals == 5000L) &&
            all(summary$markers == 50000L) && all(summary$all_gates_pass),
          "Production summary is incomplete or contains a failed panel")
qg_assert(all(acceptance$pass), "Acceptance table contains a failed maintained gate")

prefixes <- setNames(file.path(root, paste0("human_1000g_", c("eur", "afr")),
                               paste0("human_1000g_", c("eur", "afr"))),
                     c("EUR", "AFR"))
bims <- list()

for (ancestry in names(prefixes)) {
  prefix <- prefixes[[ancestry]]
  ids <- sprintf("MOSAIC_%s_%05d", ancestry, seq_len(5000L))
  check <- qg_validate_plink_files(prefix, selected, ids)
  bims[[ancestry]] <- check$bim

  bed_header <- readBin(paste0(prefix, ".bed"), what = "raw", n = 3L)
  qg_assert(identical(as.integer(bed_header), c(108L, 27L, 1L)),
            paste(ancestry, "BED is not SNP-major PLINK1"))

  recorded <- checksums[checksums$ancestry == ancestry, , drop = FALSE]
  qg_assert(nrow(recorded) == 3L, paste(ancestry, "checksum rows are incomplete"))
  for (extension in c("bed", "bim", "fam")) {
    row <- recorded[recorded$extension == extension, , drop = FALSE]
    path <- file.path(root, row$path)
    qg_assert(file.info(path)$size == row$bytes,
              paste(ancestry, extension, "size differs from manifest"))
    qg_assert(identical(qg_sha256(path), row$sha256),
              paste(ancestry, extension, "SHA-256 differs from manifest"))
  }

  focus <- unique(round(seq(1, nrow(selected), length.out = 100L)))
  work <- file.path(root, "cache", "human_validation",
                    paste0(tolower(ancestry), "-", Sys.getpid()))
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  id_path <- file.path(work, "markers.txt")
  data.table::fwrite(data.frame(ID = selected$ID[focus]), id_path,
                     col.names = FALSE)
  out <- file.path(work, "focused")
  log <- system2(plink2, c("--bfile", prefix, "--extract", id_path,
                          "--make-bed", "--out", out),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  qg_assert(is.null(status) || status == 0L,
            paste(ancestry, "focused PLINK2 extraction failed"))
  genio_back <- t(genio::read_plink(out, verbose = FALSE)$X)
  qg_assert(identical(dim(genio_back), c(5000L, 100L)) &&
              qg_valid_hard_calls(genio_back),
            paste(ancestry, "focused genio read-back failed"))

  qgg_dir <- file.path(work, "qgg")
  dir.create(qgg_dir)
  old <- setwd(qgg_dir)
  glist <- qgg::gprep(study = paste0("validate_", tolower(ancestry)),
                      bedfiles = paste0(out, ".bed"),
                      bimfiles = paste0(out, ".bim"),
                      famfiles = paste0(out, ".fam"), overwrite = TRUE)
  setwd(old)
  qgg_back <- qgg::getG(Glist = glist, chr = 1, cls = seq_len(100L),
                        impute = FALSE, scale = FALSE)
  qg_assert(identical(dim(qgg_back), dim(genio_back)) &&
              identical(rownames(qgg_back), rownames(genio_back)) &&
              identical(colnames(qgg_back), colnames(genio_back)) &&
              isTRUE(all.equal(unname(qgg_back), unname(genio_back),
                               tolerance = 0, check.attributes = FALSE)),
            paste(ancestry, "focused qgg/genio read-backs differ"))
}

qg_assert(identical(bims$EUR, bims$AFR),
          "EUR and AFR marker metadata/order are not identical")
cat("status=PASS\n")
cat("panels=EUR,AFR\n")
cat("dimensions_each=5000x50000\n")
cat("bed_bytes_each=62500003\n")
cat("checksums_match=TRUE\n")
cat("marker_metadata_identical=TRUE\n")
cat("focused_genio_qgg_roundtrip=TRUE\n")

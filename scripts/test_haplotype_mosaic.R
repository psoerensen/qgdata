#!/usr/bin/env Rscript

required <- c("digest", "genio")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "))

find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "qgdata.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find qgdata root")
    path <- parent
  }
}
root <- find_root()
source(file.path(root, "scripts", "human_1000g_common.R"), local = TRUE)
source(file.path(root, "scripts", "human_panel_validation_common.R"), local = TRUE)
source(file.path(root, "scripts", "simulate_haplotype_mosaic.R"), local = TRUE)

assert_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ },
                      error = function(e) conditionMessage(e))
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}

H <- rbind(d1 = c(0, 0, 1, 1, 0), d2 = c(1, 1, 0, 0, 1),
           d3 = c(0, 1, 0, 1, 0))
colnames(H) <- paste0("m", seq_len(ncol(H)))
map <- c(0, 0.001, 0.002, 0.004, 0.010)
a <- simulate_haplotype_mosaic(H, map, 25, 50, 12345, TRUE)
b <- simulate_haplotype_mosaic(H, map, 25, 50, 12345, TRUE)
stopifnot(identical(a, b), identical(colnames(a$haplotypes), colnames(H)),
          all(a$haplotypes %in% 0:1), identical(dim(a$haplotypes), c(25L, 5L)))

zero <- simulate_haplotype_mosaic(H, map, 50, 0, 12346, TRUE)
stopifnot(all(zero$diagnostics$donor_switches == 0L),
          all(apply(zero$haplotypes, 1L, function(x)
            any(apply(H, 1L, function(donor) all(donor == x))))))

probability <- 1 - exp(-10 * 0.1)
empirical <- simulate_haplotype_mosaic(H[, 1:2], c(0, 0.1), 5000, 10,
                                       12347, TRUE)
observed <- mean(empirical$diagnostics$donor_switches)
stopifnot(abs(observed - probability) < 0.03)
switched <- empirical$segments[empirical$segments$synthetic_haplotype %in%
  empirical$diagnostics$synthetic_haplotype[
    empirical$diagnostics$donor_switches > 0], ]
stopifnot(all(vapply(split(switched$donor, switched$synthetic_haplotype),
                     function(x) all(x[-1L] != x[-length(x)]), logical(1))))

diploid <- simulate_diploid_mosaic(H, map, 12, 25, 12348, "TEST")
stopifnot(identical(dim(diploid$genotypes), c(12L, 5L)),
          all(diploid$genotypes %in% 0:2),
          identical(colnames(diploid$genotypes), colnames(H)))

assert_error(simulate_haplotype_mosaic(matrix(c(0, 2), 1), c(0, 1), 1, 0, 1),
             "binary 0/1")
assert_error(simulate_haplotype_mosaic(H, rev(map), 1, 0, 1),
             "monotonically non-decreasing")
bad <- H
bad[1, 1] <- NA
assert_error(simulate_haplotype_mosaic(bad, map, 1, 0, 1), "complete phased")

known_segments <- data.frame(
  synthetic_haplotype = rep(sprintf("HAP%05d", 1:4), each = 2),
  donor = c(1, 2, 1, 3, 1, 4, 2, 4),
  start_marker = rep(c(1L, 4L), 4), end_marker = rep(c(3L, 5L), 4),
  stringsAsFactors = FALSE)
known_donor <- segments_to_donor_matrix(known_segments, 4L, 5L)
known_overlap <- individual_pair_donor_overlap(
  known_donor, c(0, .01, .02, .03, .04), c(100L, 200L, 300L, 400L, 500L),
  c("A", "B"))
stopifnot(nrow(known_overlap) == 1L,
          identical(known_overlap$total_same_donor_fraction, 0.3),
          identical(known_overlap$maximum_haplotype_pair_overlap, 0.6),
          identical(known_overlap$longest_same_donor_morgan, 0.02),
          identical(known_overlap$longest_same_donor_bp, 200))
known_summary <- sampled_individual_donor_overlap_summary(known_donor,
  c(0, .01, .02, .03, .04), c(100L, 200L, 300L, 400L, 500L), 1L, 44L)
stopifnot(known_summary$sampled_pairs == 1L,
          identical(known_summary$total_same_donor_mean, 0.3),
          identical(known_summary$maximum_haplotype_overlap_mean, 0.6),
          identical(known_summary$longest_tract_morgan_median, 0.02),
          identical(known_summary$longest_tract_bp_median, 200))

rate_summary <- data.frame(copy_rate = c(5, 10, 25),
  ld_pass = c(TRUE, TRUE, FALSE), relatedness_pass = c(FALSE, TRUE, TRUE),
  other_gates_pass = TRUE)
stopifnot(identical(select_highest_mosaic_rate(rate_summary), 10),
          is.na(select_highest_mosaic_rate(transform(rate_summary,
            other_gates_pass = FALSE))))

small_G <- rbind(c(0, 0, 1, 2), c(0, 1, 1, 2), c(1, 1, 2, 2),
                 c(2, 2, 1, 2), c(2, 1, 0, 2))
full_C <- suppressWarnings(stats::cor(small_G))
full_C[!is.finite(full_C)] <- 0
diag(full_C) <- 1
sample_eigen <- qg_eigen_metrics_sample_space(small_G)
stopifnot(abs(sample_eigen[["effective_rank"]] -
                ncol(full_C)^2 / sum(full_C^2)) < 1e-10,
          abs(sample_eigen[["largest_eigenvalue"]] -
                max(eigen(full_C, symmetric = TRUE, only.values = TRUE)$values)) < 1e-10)

block_seed <- 98765L
one_state <- initialize_haplotype_mosaic_state(20L, nrow(H), block_seed,
                                               map[[1L]], 100)
one <- simulate_haplotype_mosaic_block(H, map, seq(100, 500, 100), 25,
                                       one_state, TRUE)
split_state <- initialize_haplotype_mosaic_state(20L, nrow(H), block_seed,
                                                 map[[1L]], 100)
left <- simulate_haplotype_mosaic_block(H[, 1:2, drop = FALSE], map[1:2],
                                        c(100, 200), 25, split_state, TRUE)
right <- simulate_haplotype_mosaic_block(H[, 3:5, drop = FALSE], map[3:5],
                                         c(300, 400, 500), 25, left$state, TRUE)
stopifnot(identical(one$genotypes, cbind(left$genotypes, right$genotypes)),
          identical(one$donors, cbind(left$donors, right$donors)),
          identical(one$state$donor, right$state$donor),
          identical(one$state$switches, right$state$switches))

test_dir <- file.path(root, "cache", "human_r_mosaic", "focused_test")
dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(test_dir, "roundtrip")
bim <- data.frame(chr = 22, id = colnames(H), posg = map * 100,
                  pos = seq(100L, by = 100L, length.out = ncol(H)),
                  ref = "A", alt = "G", stringsAsFactors = FALSE)
fam <- qg_plink_fam(rownames(diploid$genotypes))
genio::write_plink(prefix, t(diploid$genotypes), bim = bim, fam = fam,
                   verbose = FALSE)
selected <- data.frame(CHROM = bim$chr, ID = bim$id, CM = bim$posg,
                       POS = bim$pos, REF = bim$ref, ALT = bim$alt,
                       stringsAsFactors = FALSE)
qg_validate_plink_files(prefix, selected, rownames(diploid$genotypes))
back <- t(genio::read_plink(prefix, verbose = FALSE)$X)
stopifnot(identical(dim(back), dim(diploid$genotypes)),
          isTRUE(all.equal(unname(back), unname(diploid$genotypes),
                           tolerance = 0, check.attributes = FALSE)))

cat("status=PASS\n")
cat("reproducibility=TRUE\n")
cat("zero_switch_copies_one_donor=TRUE\n")
cat("expected_switch_probability=", probability, "\n", sep = "")
cat("empirical_switch_probability=", observed, "\n", sep = "")
cat("invalid_inputs_rejected=TRUE\n")
cat("donor_overlap_calculation=TRUE\n")
cat("sampled_donor_summary=TRUE\n")
cat("deterministic_highest_rate_selection=TRUE\n")
cat("sample_space_eigenvalues_exact=TRUE\n")
cat("blockwise_generation_exact=TRUE\n")
cat("plink_roundtrip_exact=TRUE\n")

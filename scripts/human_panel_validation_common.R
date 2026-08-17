# Shared metrics for sample-size-matched human genotype-panel validation.

qg_validation_distance_breaks <- c(-Inf, 1e4, 5e4, 1e5, 5e5, 1e6, Inf)
qg_validation_distance_labels <- c("0-10kb", "10-50kb", "50-100kb",
                                   "100-500kb", "500kb-1Mb", ">1Mb")
qg_validation_distance_names <- c("0_10kb", "10_50kb", "50_100kb",
                                  "100_500kb", "500kb_1Mb", "over_1Mb")

qg_ld_by_distance <- function(G, positions) {
  qg_assert(is.matrix(G) && ncol(G) == length(positions),
            "LD matrix and position count differ")
  C <- suppressWarnings(stats::cor(G))
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  upper <- which(upper.tri(C), arr.ind = TRUE)
  distances <- abs(positions[upper[, 1]] - positions[upper[, 2]])
  bins <- cut(distances, breaks = qg_validation_distance_breaks,
              labels = qg_validation_distance_labels, right = FALSE)
  r2 <- C[upper]^2
  do.call(rbind, lapply(seq_along(qg_validation_distance_labels), function(k) {
    take <- which(bins == qg_validation_distance_labels[[k]])
    data.frame(distance_bin = qg_validation_distance_labels[[k]],
      distance_name = qg_validation_distance_names[[k]],
      marker_pairs = length(take),
      mean_r2 = if (length(take)) mean(r2[take]) else NA_real_,
      median_r2 = if (length(take)) stats::median(r2[take]) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

qg_row_hashes <- function(G) vapply(seq_len(nrow(G)), function(i)
  digest::digest(G[i, ], algo = "sha256", serialize = TRUE), character(1))

qg_column_hashes <- function(G) vapply(seq_len(ncol(G)), function(j)
  digest::digest(G[, j], algo = "sha256", serialize = TRUE), character(1))

qg_identical_marker_pairs <- function(G) {
  counts <- table(qg_column_hashes(G))
  as.double(sum(counts * (counts - 1) / 2))
}

qg_largest_eigenvalue <- function(C, iterations = 150L) {
  v <- rep(1 / sqrt(ncol(C)), ncol(C))
  for (i in seq_len(iterations)) {
    next_v <- as.vector(C %*% v)
    norm <- sqrt(sum(next_v^2))
    if (!is.finite(norm) || norm == 0) return(0)
    v <- next_v / norm
  }
  as.numeric(crossprod(v, C %*% v))
}

# KING-robust kinship estimator: (HetHet - 2*IBS0)/(Het_i + Het_j).
qg_king_robust_kinship <- function(G) {
  heterozygous <- G == 1
  hom_ref <- G == 0
  hom_alt <- G == 2
  both_het <- tcrossprod(heterozygous)
  ibs0 <- tcrossprod(hom_ref, hom_alt) + tcrossprod(hom_alt, hom_ref)
  denominator <- outer(rowSums(heterozygous), rowSums(heterozygous), "+")
  kinship <- (both_het - 2 * ibs0) / denominator
  diag(kinship) <- NA_real_
  kinship
}

qg_kinship_metrics <- function(G) {
  K <- qg_king_robust_kinship(G)
  values <- K[upper.tri(K)]
  values <- values[is.finite(values)]
  c(
    kinship_pairs = length(values),
    kinship_median = stats::median(values),
    kinship_q95 = unname(stats::quantile(values, 0.95)),
    kinship_max = max(values),
    duplicate_MZ_pairs = sum(values >= 2^(-1.5)),
    first_degree_pairs = sum(values >= 2^(-2.5) & values < 2^(-1.5)),
    second_degree_pairs = sum(values >= 2^(-3.5) & values < 2^(-2.5))
  )
}

qg_kinship_distribution <- function(G) {
  K <- qg_king_robust_kinship(G)
  values <- K[upper.tri(K)]
  values <- values[is.finite(values)]
  probabilities <- c(0, .025, .25, .5, .75, .9, .95, .975, .99, .995, .999, 1)
  quantiles <- stats::quantile(values, probabilities, names = FALSE)
  names(quantiles) <- c("minimum", "q025", "q250", "median", "q750", "q900",
                        "q950", "q975", "q990", "q995", "q999", "maximum")
  c(pairs = length(values), quantiles,
    duplicate_MZ_pairs = sum(values >= 2^(-1.5)),
    first_degree_pairs = sum(values >= 2^(-2.5) & values < 2^(-1.5)),
    second_degree_pairs = sum(values >= 2^(-3.5) & values < 2^(-2.5)))
}

qg_eigen_metrics_sample_space <- function(G) {
  centered <- sweep(G, 2L, colMeans(G), "-")
  norm <- sqrt(colSums(centered^2))
  polymorphic <- is.finite(norm) & norm > 0
  Z <- centered[, polymorphic, drop = FALSE]
  Z <- sweep(Z, 2L, norm[polymorphic], "/")
  eigenvalues <- eigen(tcrossprod(Z), symmetric = TRUE, only.values = TRUE)$values
  eigenvalues[eigenvalues < 0 & eigenvalues > -1e-8] <- 0
  monomorphic <- sum(!polymorphic)
  sum_squares <- sum(eigenvalues^2) + monomorphic
  c(effective_rank = ncol(G)^2 / sum_squares,
    largest_eigenvalue = max(c(eigenvalues, if (monomorphic) 1 else 0)))
}

qg_local_ld_metrics <- function(G, positions, max_distance = 1e6) {
  centered <- sweep(G, 2L, colMeans(G), "-")
  norm <- sqrt(colSums(centered^2))
  Z <- matrix(0, nrow(G), ncol(G))
  keep <- is.finite(norm) & norm > 0
  Z[, keep] <- sweep(centered[, keep, drop = FALSE], 2L, norm[keep], "/")
  ld_scores <- numeric(ncol(G))
  r2_by_bin <- vector("list", length(qg_validation_distance_labels))
  for (j in seq_len(ncol(G) - 1L)) {
    end <- findInterval(positions[[j]] + max_distance, positions)
    if (end <= j) next
    index <- (j + 1L):end
    correlation <- as.numeric(crossprod(Z[, j], Z[, index, drop = FALSE]))
    r2 <- correlation^2
    ld_scores[[j]] <- ld_scores[[j]] + sum(r2)
    ld_scores[index] <- ld_scores[index] + r2
    bins <- cut(positions[index] - positions[[j]],
      breaks = qg_validation_distance_breaks,
      labels = qg_validation_distance_labels, right = FALSE)
    for (k in seq_along(r2_by_bin)) {
      take <- bins == qg_validation_distance_labels[[k]]
      if (any(take)) r2_by_bin[[k]][[length(r2_by_bin[[k]]) + 1L]] <- r2[take]
    }
  }
  rows <- do.call(rbind, lapply(seq_along(qg_validation_distance_labels), function(k) {
    values <- if (k <= length(r2_by_bin)) unlist(r2_by_bin[[k]], use.names = FALSE)
              else numeric()
    data.frame(distance_bin = qg_validation_distance_labels[[k]],
      distance_name = qg_validation_distance_names[[k]], marker_pairs = length(values),
      mean_r2 = if (length(values)) mean(values) else NA_real_,
      median_r2 = if (length(values)) stats::median(values) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  list(ld_by_distance = rows, median_ld_score = stats::median(ld_scores))
}

qg_dataset_metrics_scalable <- function(G, positions, full_reference_maf,
                                        pca_reference = NULL,
                                        max_ld_distance = 1e6) {
  qg_assert(qg_valid_hard_calls(G),
            "Validation matrix is not complete diploid 0/1/2 hard calls")
  qg_assert(length(positions) == ncol(G), "Position count differs from genotype markers")
  qg_assert(length(full_reference_maf) == ncol(G),
            "Reference-MAF count differs from genotype markers")
  maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
  local <- qg_local_ld_metrics(G, positions, max_ld_distance)
  eigen_metrics <- qg_eigen_metrics_sample_space(G)
  values <- c(
    maf_correlation_with_full_reference = stats::cor(maf, full_reference_maf),
    mean_maf = mean(maf), median_maf = stats::median(maf),
    mean_maf_absolute_difference = mean(abs(maf - full_reference_maf)),
    median_maf_absolute_difference = stats::median(abs(maf - full_reference_maf)),
    max_maf_absolute_difference = max(abs(maf - full_reference_maf)),
    genotype_missingness = mean(is.na(G)), heterozygosity = mean(G == 1),
    median_ld_score = local$median_ld_score,
    eigen_metrics,
    monomorphic_markers = sum(vapply(seq_len(ncol(G)), function(j)
      length(unique(G[, j])) < 2L, logical(1))),
    duplicated_individuals = sum(duplicated(qg_row_hashes(G))),
    identical_genotype_marker_pairs = qg_identical_marker_pairs(G),
    qg_kinship_metrics(G)
  )
  if (!is.null(pca_reference)) values <- c(values, qg_pca_metrics(G, pca_reference))
  for (k in seq_along(qg_validation_distance_labels)) {
    values[[paste0("pairs_", qg_validation_distance_names[[k]])]] <-
      local$ld_by_distance$marker_pairs[[k]]
    values[[paste0("mean_r2_", qg_validation_distance_names[[k]])]] <-
      local$ld_by_distance$mean_r2[[k]]
    values[[paste0("median_r2_", qg_validation_distance_names[[k]])]] <-
      local$ld_by_distance$median_r2[[k]]
  }
  values
}

qg_prepare_pca_reference <- function(R, components = 10L) {
  variances <- apply(R, 2L, stats::var)
  keep <- is.finite(variances) & variances > 0
  qg_assert(sum(keep) >= 2L, "Too few polymorphic reference markers for PCA")
  rank <- min(as.integer(components), nrow(R) - 1L, sum(keep))
  fit <- stats::prcomp(R[, keep, drop = FALSE], center = TRUE, scale. = TRUE,
                       rank. = rank)
  list(keep = keep, center = fit$center, scale = fit$scale,
       rotation = fit$rotation[, seq_len(rank), drop = FALSE],
       reference_scores = fit$x[, seq_len(rank), drop = FALSE])
}

qg_pca_metrics <- function(G, pca_reference) {
  X <- sweep(G[, pca_reference$keep, drop = FALSE], 2L,
             pca_reference$center, "-")
  X <- sweep(X, 2L, pca_reference$scale, "/")
  scores <- X %*% pca_reference$rotation
  reference_scores <- pca_reference$reference_scores
  reference_sd <- apply(reference_scores, 2L, stats::sd)
  centroid_z <- (colMeans(scores) - colMeans(reference_scores)) / reference_sd
  variance_ratio <- apply(scores, 2L, stats::var) /
    apply(reference_scores, 2L, stats::var)
  c(
    pca_components = ncol(scores),
    pca_centroid_rms_z = sqrt(mean(centroid_z^2)),
    pca_centroid_max_abs_z = max(abs(centroid_z)),
    pca_variance_ratio_median = stats::median(variance_ratio)
  )
}

qg_dataset_metrics <- function(G, positions, full_reference_maf,
                               pca_reference = NULL) {
  qg_assert(qg_valid_hard_calls(G),
            "Validation matrix is not complete diploid 0/1/2 hard calls")
  qg_assert(length(positions) == ncol(G), "Position count differs from genotype markers")
  qg_assert(length(full_reference_maf) == ncol(G),
            "Reference-MAF count differs from genotype markers")
  maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
  maf_abs_difference <- abs(maf - full_reference_maf)
  C <- suppressWarnings(stats::cor(G))
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  ld_by_distance <- qg_ld_by_distance(G, positions)
  within_1mb <- abs(outer(positions, positions, "-")) <= 1e6
  diag(within_1mb) <- FALSE
  ld_scores <- rowSums(C^2 * within_1mb)
  values <- c(
    maf_correlation_with_full_reference = stats::cor(maf, full_reference_maf),
    mean_maf = mean(maf), median_maf = stats::median(maf),
    mean_maf_absolute_difference = mean(maf_abs_difference),
    median_maf_absolute_difference = stats::median(maf_abs_difference),
    max_maf_absolute_difference = max(maf_abs_difference),
    genotype_missingness = mean(is.na(G)), heterozygosity = mean(G == 1),
    median_ld_score = stats::median(ld_scores),
    effective_rank = ncol(C)^2 / sum(C^2),
    largest_eigenvalue = qg_largest_eigenvalue(C),
    monomorphic_markers = sum(vapply(seq_len(ncol(G)), function(j)
      length(unique(G[, j])) < 2L, logical(1))),
    duplicated_individuals = sum(duplicated(qg_row_hashes(G))),
    identical_genotype_marker_pairs = qg_identical_marker_pairs(G),
    qg_kinship_metrics(G)
  )
  if (!is.null(pca_reference)) values <- c(values, qg_pca_metrics(G, pca_reference))
  for (k in seq_along(qg_validation_distance_labels)) {
    values[[paste0("pairs_", qg_validation_distance_names[[k]])]] <-
      ld_by_distance$marker_pairs[[k]]
    values[[paste0("mean_r2_", qg_validation_distance_names[[k]])]] <-
      ld_by_distance$mean_r2[[k]]
    values[[paste0("median_r2_", qg_validation_distance_names[[k]])]] <-
      ld_by_distance$median_r2[[k]]
  }
  values
}

qg_matched_reference_indices <- function(n_reference, n_subsamples = 20L,
                                         n_match = 200L, seed) {
  qg_assert(n_reference >= n_match, "Reference panel is smaller than matched sample size")
  set.seed(as.integer(seed))
  replicate(n_subsamples, sort(sample.int(n_reference, n_match)), simplify = FALSE)
}

qg_metric_distribution_row <- function(candidate, reference_values,
                                       metadata) {
  data.frame(
    metadata,
    metric = names(candidate),
    synthetic_value = unname(candidate),
    reference_mean = colMeans(reference_values, na.rm = TRUE),
    reference_sd = apply(reference_values, 2L, stats::sd, na.rm = TRUE),
    reference_median = apply(reference_values, 2L, stats::median, na.rm = TRUE),
    reference_q025 = apply(reference_values, 2L, stats::quantile, 0.025,
                           na.rm = TRUE),
    reference_q975 = apply(reference_values, 2L, stats::quantile, 0.975,
                           na.rm = TRUE),
    synthetic_empirical_percentile = vapply(seq_along(candidate), function(i)
      mean(reference_values[, i] <= candidate[[i]], na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

qg_gate_specification <- function() {
  data.frame(
    gate = c("dimensions_hard_calls", "marker_integrity", "duplicate_individuals",
             "maf_correlation", "mean_maf_difference", "heterozygosity_difference",
             "ld_below_100kb", "effective_rank", "largest_eigenvalue",
             "close_related_pairs", "pca_displacement", "replicate_consistency"),
    pass_rule = c(
      "exact requested dimensions; complete diploid 0/1/2 calls",
      "exact unique marker IDs/order/positions/alleles; no omissions",
      "0 duplicated synthetic individuals", ">=0.98",
      "absolute synthetic minus reference mean <=0.01",
      "absolute synthetic minus reference mean <=0.01",
      "for each <100kb bin with >=10 pairs: replicate median within 20% of reference mean and preferably within reference 95% interval",
      "replicate median within 20% of matched-reference mean",
      "replicate median within 20% of matched-reference mean",
      "0 duplicate/MZ and 0 first-degree synthetic pairs",
      "median PCA centroid RMS z <=0.25 pass; >0.25 to 0.5 borderline; >0.5 fail",
      ">=5 deterministic replicates and no gate depends on one favourable seed"
    ), stringsAsFactors = FALSE)
}

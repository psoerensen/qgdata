# Transparent reference-conditioned haplotype-mosaic simulation.

validate_haplotype_mosaic_inputs <- function(H, map_morgan, n, copy_rate, seed) {
  if (!is.matrix(H) || !is.numeric(H) && !is.integer(H) && !is.logical(H))
    stop("H must be a numeric, integer, or logical matrix", call. = FALSE)
  if (nrow(H) < 1L || ncol(H) < 1L)
    stop("H must contain at least one donor and one marker", call. = FALSE)
  if (anyNA(H) || any(!is.finite(H)) || any(!H %in% 0:1))
    stop("H must contain only complete phased binary 0/1 alleles", call. = FALSE)
  if (!is.numeric(map_morgan) || length(map_morgan) != ncol(H) ||
      anyNA(map_morgan) || any(!is.finite(map_morgan)))
    stop("map_morgan must be finite and have one value per marker", call. = FALSE)
  if (any(diff(map_morgan) < 0))
    stop("map_morgan must be monotonically non-decreasing", call. = FALSE)
  if (length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n))
    stop("n must be a positive integer", call. = FALSE)
  if (length(copy_rate) != 1L || is.na(copy_rate) || !is.finite(copy_rate) ||
      copy_rate < 0)
    stop("copy_rate must be one finite non-negative value", call. = FALSE)
  if (length(seed) != 1L || is.na(seed) || seed != as.integer(seed))
    stop("seed must be one valid R integer", call. = FALSE)
  if (nrow(H) < 2L && copy_rate > 0)
    stop("At least two donors are required when copy_rate is positive", call. = FALSE)
  invisible(TRUE)
}

simulate_haplotype_mosaic <- function(H, map_morgan, n, copy_rate, seed,
                                      return_segments = FALSE) {
  validate_haplotype_mosaic_inputs(H, map_morgan, n, copy_rate, seed)
  H <- matrix(as.integer(H), nrow = nrow(H), ncol = ncol(H),
              dimnames = dimnames(H))
  n <- as.integer(n)
  seed <- as.integer(seed)
  delta <- c(0, diff(map_morgan))
  switch_probability <- -expm1(-copy_rate * delta)
  switch_probability[1L] <- 0
  set.seed(seed)

  result <- matrix(0L, nrow = n, ncol = ncol(H),
                   dimnames = list(sprintf("HAP%05d", seq_len(n)), colnames(H)))
  switches <- integer(n)
  distinct_donors <- integer(n)
  segments <- vector("list", n)

  for (i in seq_len(n)) {
    donor <- sample.int(nrow(H), 1L)
    donors_used <- donor
    start <- 1L
    hap_segments <- list()
    result[i, 1L] <- H[donor, 1L]
    if (ncol(H) > 1L) for (j in 2:ncol(H)) {
      if (stats::runif(1L) < switch_probability[[j]]) {
        hap_segments[[length(hap_segments) + 1L]] <- c(
          donor = donor, start_marker = start, end_marker = j - 1L)
        next_donor <- sample.int(nrow(H) - 1L, 1L)
        if (next_donor >= donor) next_donor <- next_donor + 1L
        donor <- next_donor
        donors_used <- c(donors_used, donor)
        switches[[i]] <- switches[[i]] + 1L
        start <- j
      }
      result[i, j] <- H[donor, j]
    }
    hap_segments[[length(hap_segments) + 1L]] <- c(
      donor = donor, start_marker = start, end_marker = ncol(H))
    distinct_donors[[i]] <- length(unique(donors_used))
    segment_matrix <- do.call(rbind, hap_segments)
    segments[[i]] <- data.frame(
      synthetic_haplotype = rownames(result)[[i]],
      segment = seq_len(nrow(segment_matrix)),
      donor = as.integer(segment_matrix[, "donor"]),
      start_marker = as.integer(segment_matrix[, "start_marker"]),
      end_marker = as.integer(segment_matrix[, "end_marker"]),
      markers = as.integer(segment_matrix[, "end_marker"] -
                             segment_matrix[, "start_marker"] + 1L),
      length_morgan = map_morgan[segment_matrix[, "end_marker"]] -
        map_morgan[segment_matrix[, "start_marker"]],
      stringsAsFactors = FALSE)
  }

  segment_table <- do.call(rbind, segments)
  rownames(segment_table) <- NULL
  diagnostics <- data.frame(
    synthetic_haplotype = rownames(result), donor_switches = switches,
    segments = switches + 1L, distinct_donors = distinct_donors,
    stringsAsFactors = FALSE)
  summary <- c(
    haplotypes = n, donors_available = nrow(H), markers = ncol(H),
    copy_rate = copy_rate, map_length_morgan = max(map_morgan) - min(map_morgan),
    switches_mean = mean(switches), switches_median = stats::median(switches),
    switches_q025 = unname(stats::quantile(switches, 0.025)),
    switches_q975 = unname(stats::quantile(switches, 0.975)),
    distinct_donors_mean = mean(distinct_donors),
    distinct_donors_median = stats::median(distinct_donors),
    segment_markers_q025 = unname(stats::quantile(segment_table$markers, 0.025)),
    segment_markers_median = stats::median(segment_table$markers),
    segment_markers_q975 = unname(stats::quantile(segment_table$markers, 0.975)),
    segment_morgan_q025 = unname(stats::quantile(segment_table$length_morgan, 0.025)),
    segment_morgan_median = stats::median(segment_table$length_morgan),
    segment_morgan_q975 = unname(stats::quantile(segment_table$length_morgan, 0.975))
  )
  attr(result, "diagnostics") <- diagnostics
  attr(result, "segment_summary") <- summary
  if (return_segments)
    return(list(haplotypes = result, diagnostics = diagnostics,
                segment_summary = summary, segments = segment_table,
                switch_probability = switch_probability))
  result
}

simulate_diploid_mosaic <- function(H, map_morgan, n, copy_rate, seed,
                                    individual_prefix = "MOSAIC") {
  simulated <- simulate_haplotype_mosaic(H, map_morgan, 2L * as.integer(n),
                                         copy_rate, seed,
                                         return_segments = TRUE)
  first <- seq.int(1L, 2L * n, by = 2L)
  G <- simulated$haplotypes[first, , drop = FALSE] +
    simulated$haplotypes[first + 1L, , drop = FALSE]
  storage.mode(G) <- "integer"
  rownames(G) <- sprintf("%s%05d", individual_prefix, seq_len(n))
  list(genotypes = G, haplotype_diagnostics = simulated$diagnostics,
       segment_summary = simulated$segment_summary,
       segments = simulated$segments)
}

segments_to_donor_matrix <- function(segments, n_haplotypes, n_markers) {
  required <- c("synthetic_haplotype", "donor", "start_marker", "end_marker")
  if (!is.data.frame(segments) || !all(required %in% names(segments)))
    stop("segments lacks required donor-provenance columns", call. = FALSE)
  if (n_haplotypes < 1L || n_markers < 1L)
    stop("n_haplotypes and n_markers must be positive", call. = FALSE)
  expected_ids <- sprintf("HAP%05d", seq_len(n_haplotypes))
  if (!setequal(unique(segments$synthetic_haplotype), expected_ids))
    stop("segments does not contain exactly the expected haplotypes", call. = FALSE)
  donor <- matrix(NA_integer_, nrow = n_haplotypes, ncol = n_markers,
                  dimnames = list(expected_ids, NULL))
  for (i in seq_len(nrow(segments))) {
    row <- segments[i, ]
    hap <- match(row$synthetic_haplotype, expected_ids)
    start <- as.integer(row$start_marker)
    end <- as.integer(row$end_marker)
    if (is.na(start) || is.na(end) || start < 1L || end > n_markers || start > end)
      stop("segments contains an invalid marker interval", call. = FALSE)
    if (any(!is.na(donor[hap, start:end])))
      stop("segments contains overlapping intervals", call. = FALSE)
    donor[hap, start:end] <- as.integer(row$donor)
  }
  if (anyNA(donor)) stop("segments does not cover every haplotype marker", call. = FALSE)
  donor
}

donor_overlap_matrices <- function(donor, map_morgan, positions,
                                   calculate_longest = TRUE) {
  if (!is.matrix(donor) || nrow(donor) %% 2L != 0L || anyNA(donor))
    stop("donor must be a complete even-row haplotype-by-marker matrix", call. = FALSE)
  if (length(map_morgan) != ncol(donor) || length(positions) != ncol(donor) ||
      any(!is.finite(map_morgan)) || any(diff(map_morgan) < 0) ||
      any(!is.finite(positions)) || any(diff(positions) < 0))
    stop("map and physical positions must be finite, ordered, and marker-aligned",
         call. = FALSE)
  nh <- nrow(donor)
  overlap_count <- matrix(0L, nh, nh)
  if (calculate_longest) {
    run_start <- matrix(0L, nh, nh)
    longest_morgan <- matrix(0, nh, nh)
    longest_bp <- matrix(0, nh, nh)
  }
  close_runs <- function(ending, end_marker) {
    if (!any(ending)) return(invisible(NULL))
    starts <- run_start[ending]
    length_morgan <- map_morgan[[end_marker]] - map_morgan[starts]
    length_bp <- positions[[end_marker]] - positions[starts]
    current_m <- longest_morgan[ending]
    current_bp <- longest_bp[ending]
    longest_morgan[ending] <<- pmax(current_m, length_morgan)
    longest_bp[ending] <<- pmax(current_bp, length_bp)
    run_start[ending] <<- 0L
    invisible(NULL)
  }
  for (j in seq_len(ncol(donor))) {
    equal <- outer(donor[, j], donor[, j], `==`)
    overlap_count <- overlap_count + equal
    if (calculate_longest) {
      starting <- equal & run_start == 0L
      run_start[starting] <- j
      ending <- !equal & run_start > 0L
      close_runs(ending, j - 1L)
    }
  }
  if (calculate_longest) close_runs(run_start > 0L, ncol(donor))
  list(overlap_fraction = overlap_count / ncol(donor),
       longest_morgan = if (calculate_longest) longest_morgan else NULL,
       longest_bp = if (calculate_longest) longest_bp else NULL)
}

individual_pair_donor_overlap <- function(donor, map_morgan, positions,
                                          individual_ids = NULL,
                                          kinship = NULL,
                                          genotype_correlation = NULL,
                                          calculate_longest = TRUE) {
  matrices <- donor_overlap_matrices(donor, map_morgan, positions,
                                     calculate_longest)
  n <- nrow(donor) / 2L
  if (is.null(individual_ids)) individual_ids <- sprintf("IND%05d", seq_len(n))
  if (length(individual_ids) != n)
    stop("individual_ids length differs from paired haplotypes", call. = FALSE)
  pairs <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
  total <- maximum <- longest_m <- longest_bp <- numeric(nrow(pairs))
  kin <- correlation <- rep(NA_real_, nrow(pairs))
  for (p in seq_len(nrow(pairs))) {
    i <- pairs[p, 1L]
    j <- pairs[p, 2L]
    hi <- c(2L * i - 1L, 2L * i)
    hj <- c(2L * j - 1L, 2L * j)
    overlap <- as.vector(matrices$overlap_fraction[hi, hj, drop = FALSE])
    total[[p]] <- mean(overlap)
    maximum[[p]] <- max(overlap)
    if (calculate_longest) {
      longest_m[[p]] <- max(matrices$longest_morgan[hi, hj, drop = FALSE])
      longest_bp[[p]] <- max(matrices$longest_bp[hi, hj, drop = FALSE])
    } else {
      longest_m[[p]] <- longest_bp[[p]] <- NA_real_
    }
    if (!is.null(kinship)) kin[[p]] <- kinship[i, j]
    if (!is.null(genotype_correlation)) correlation[[p]] <- genotype_correlation[i, j]
  }
  data.frame(individual_1 = individual_ids[pairs[, 1L]],
    individual_2 = individual_ids[pairs[, 2L]],
    total_same_donor_fraction = total,
    maximum_haplotype_pair_overlap = maximum,
    longest_same_donor_morgan = longest_m,
    longest_same_donor_bp = longest_bp,
    kinship = kin, genotype_correlation = correlation,
    stringsAsFactors = FALSE)
}

select_highest_mosaic_rate <- function(rate_summary) {
  required <- c("copy_rate", "ld_pass", "relatedness_pass", "other_gates_pass")
  if (!is.data.frame(rate_summary) || !all(required %in% names(rate_summary)))
    stop("rate_summary lacks required selection columns", call. = FALSE)
  eligible <- with(rate_summary, is.finite(copy_rate) & ld_pass &
                     relatedness_pass & other_gates_pass)
  if (!any(eligible)) return(NA_real_)
  max(rate_summary$copy_rate[eligible])
}

# Marker-major implementation of the same copying process.  The state and the
# global R RNG stream are carried across calls, so large panels can be written
# in marker blocks without materializing the complete genotype matrix.
initialize_haplotype_mosaic_state <- function(n_haplotypes, n_donors, seed,
                                               first_morgan, first_bp) {
  if (length(n_haplotypes) != 1L || n_haplotypes < 2L ||
      n_haplotypes != as.integer(n_haplotypes))
    stop("n_haplotypes must be an integer >=2", call. = FALSE)
  if (length(n_donors) != 1L || n_donors < 2L ||
      n_donors != as.integer(n_donors))
    stop("n_donors must be an integer >=2", call. = FALSE)
  if (!is.finite(first_morgan) || !is.finite(first_bp))
    stop("initial map positions must be finite", call. = FALSE)
  set.seed(as.integer(seed))
  list(donor = sample.int(as.integer(n_donors), as.integer(n_haplotypes),
                          replace = TRUE),
       n_donors = as.integer(n_donors), switches = integer(n_haplotypes),
       segment_start_morgan = rep(first_morgan, n_haplotypes),
       segment_start_bp = rep(first_bp, n_haplotypes),
       segment_morgan = numeric(), segment_bp = numeric(),
       previous_morgan = NA_real_, previous_bp = NA_real_,
       markers_generated = 0L)
}

simulate_haplotype_mosaic_block <- function(H, map_morgan, positions,
                                             copy_rate, state,
                                             return_donors = FALSE) {
  if (!is.matrix(H) || anyNA(H) || any(!H %in% 0:1))
    stop("H must be a complete phased binary 0/1 matrix", call. = FALSE)
  if (length(map_morgan) != ncol(H) || length(positions) != ncol(H) ||
      any(!is.finite(map_morgan)) || any(!is.finite(positions)) ||
      any(diff(map_morgan) < 0) || any(diff(positions) < 0))
    stop("block map and positions must be finite, ordered, and marker-aligned",
         call. = FALSE)
  if (!is.list(state) || length(state$donor) < 2L ||
      state$n_donors != nrow(H) || any(state$donor < 1L | state$donor > nrow(H)))
    stop("mosaic state is incompatible with H", call. = FALSE)
  if (length(copy_rate) != 1L || !is.finite(copy_rate) || copy_rate < 0)
    stop("copy_rate must be finite and non-negative", call. = FALSE)
  if (state$markers_generated > 0L &&
      (map_morgan[[1L]] < state$previous_morgan || positions[[1L]] < state$previous_bp))
    stop("block starts before the preceding marker", call. = FALSE)

  n_haplotypes <- length(state$donor)
  haplotypes <- matrix(0L, n_haplotypes, ncol(H),
    dimnames = list(sprintf("HAP%05d", seq_len(n_haplotypes)), colnames(H)))
  donor_matrix <- if (return_donors) matrix(0L, n_haplotypes, ncol(H)) else NULL
  for (j in seq_len(ncol(H))) {
    delta <- if (j > 1L) map_morgan[[j]] - map_morgan[[j - 1L]] else if (
      state$markers_generated > 0L) map_morgan[[j]] - state$previous_morgan else 0
    switch <- if (state$markers_generated == 0L && j == 1L) {
      rep(FALSE, n_haplotypes)
    } else {
      stats::runif(n_haplotypes) < -expm1(-copy_rate * delta)
    }
    if (any(switch)) {
      state$segment_morgan <- c(state$segment_morgan,
        map_morgan[[j]] - state$segment_start_morgan[switch])
      state$segment_bp <- c(state$segment_bp,
        positions[[j]] - state$segment_start_bp[switch])
      current <- state$donor[switch]
      next_donor <- sample.int(state$n_donors - 1L, sum(switch), replace = TRUE)
      next_donor[next_donor >= current] <- next_donor[next_donor >= current] + 1L
      state$donor[switch] <- next_donor
      state$switches[switch] <- state$switches[switch] + 1L
      state$segment_start_morgan[switch] <- map_morgan[[j]]
      state$segment_start_bp[switch] <- positions[[j]]
    }
    haplotypes[, j] <- H[cbind(state$donor, rep.int(j, n_haplotypes))]
    if (return_donors) donor_matrix[, j] <- state$donor
    state$previous_morgan <- map_morgan[[j]]
    state$previous_bp <- positions[[j]]
    state$markers_generated <- state$markers_generated + 1L
  }
  first <- seq.int(1L, n_haplotypes, by = 2L)
  genotypes <- haplotypes[first, , drop = FALSE] +
    haplotypes[first + 1L, , drop = FALSE]
  storage.mode(genotypes) <- "integer"
  list(genotypes = genotypes, donors = donor_matrix, state = state)
}

finalize_haplotype_mosaic_state <- function(state) {
  state$segment_morgan <- c(state$segment_morgan,
    state$previous_morgan - state$segment_start_morgan)
  state$segment_bp <- c(state$segment_bp,
    state$previous_bp - state$segment_start_bp)
  c(haplotypes = length(state$donor), markers = state$markers_generated,
    switches_mean = mean(state$switches),
    switches_median = stats::median(state$switches),
    switches_q025 = unname(stats::quantile(state$switches, .025)),
    switches_q975 = unname(stats::quantile(state$switches, .975)),
    segment_morgan_q025 = unname(stats::quantile(state$segment_morgan, .025)),
    segment_morgan_median = stats::median(state$segment_morgan),
    segment_morgan_q975 = unname(stats::quantile(state$segment_morgan, .975)),
    segment_bp_q025 = unname(stats::quantile(state$segment_bp, .025)),
    segment_bp_median = stats::median(state$segment_bp),
    segment_bp_q975 = unname(stats::quantile(state$segment_bp, .975)))
}

sampled_individual_donor_overlap_summary <- function(donor, map_morgan, positions,
                                                       n_pairs = 2000L, seed) {
  if (!is.matrix(donor) || nrow(donor) %% 2L != 0L || anyNA(donor) ||
      length(map_morgan) != ncol(donor) || length(positions) != ncol(donor))
    stop("donor matrix/map inputs are incompatible", call. = FALSE)
  n_individuals <- nrow(donor) / 2L
  if (n_pairs < 1L || n_pairs > n_individuals * (n_individuals - 1L) / 2L)
    stop("n_pairs is outside the available individual-pair range", call. = FALSE)
  set.seed(as.integer(seed))
  pair_key <- character()
  while (length(pair_key) < n_pairs) {
    sampled_i <- sample.int(n_individuals, max(1000L, n_pairs), replace = TRUE)
    sampled_j <- sample.int(n_individuals, max(1000L, n_pairs), replace = TRUE)
    keep <- sampled_i != sampled_j
    key <- paste(pmin(sampled_i[keep], sampled_j[keep]),
                 pmax(sampled_i[keep], sampled_j[keep]), sep = ":")
    pair_key <- unique(c(pair_key, key))
    pair_key <- pair_key[seq_len(min(n_pairs, length(pair_key)))]
  }
  pieces <- do.call(rbind, strsplit(pair_key, ":", fixed = TRUE))
  pair_i <- as.integer(pieces[, 1L]); pair_j <- as.integer(pieces[, 2L])
  longest <- function(equal) {
    runs <- rle(equal); ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L; take <- which(runs$values)
    if (!length(take)) return(c(morgan = 0, bp = 0))
    lengths <- map_morgan[ends[take]] - map_morgan[starts[take]]
    best <- take[[which.max(lengths)]]
    c(morgan = max(lengths), bp = positions[ends[best]] - positions[starts[best]])
  }
  total <- maximum <- longest_m <- longest_bp <- numeric(n_pairs)
  for (p in seq_len(n_pairs)) {
    hi <- c(2L * pair_i[[p]] - 1L, 2L * pair_i[[p]])
    hj <- c(2L * pair_j[[p]] - 1L, 2L * pair_j[[p]])
    overlap <- numeric(4L); lengths <- vector("list", 4L); cursor <- 1L
    for (a in hi) for (b in hj) {
      equal <- donor[a, ] == donor[b, ]
      overlap[[cursor]] <- mean(equal); lengths[[cursor]] <- longest(equal)
      cursor <- cursor + 1L
    }
    total[[p]] <- mean(overlap); maximum[[p]] <- max(overlap)
    longest_m[[p]] <- max(vapply(lengths, `[[`, numeric(1), "morgan"))
    longest_bp[[p]] <- max(vapply(lengths, `[[`, numeric(1), "bp"))
  }
  data.frame(sampled_pairs = n_pairs,
    total_same_donor_mean = mean(total),
    total_same_donor_q995 = unname(stats::quantile(total, .995)),
    maximum_haplotype_overlap_mean = mean(maximum),
    maximum_haplotype_overlap_q995 = unname(stats::quantile(maximum, .995)),
    longest_tract_morgan_median = stats::median(longest_m),
    longest_tract_morgan_q995 = unname(stats::quantile(longest_m, .995)),
    longest_tract_bp_median = stats::median(longest_bp),
    longest_tract_bp_q995 = unname(stats::quantile(longest_bp, .995)),
    stringsAsFactors = FALSE)
}

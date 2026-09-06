# Shared analytical functions for empirical and simulation benchmarks.

mx_log("Stage 02: loading analytical helpers.")
mx_require_packages(c("igraph", "Matrix", "RSpectra", "ggplot2", "reshape2"))
mx_assert_files(MX_PATHS$graphon_reference, "Graphon reference implementation")

# The reviewed graphon estimator is kept as a separate source file so that the
# MethodsX implementation can be compared line by line with the related study.
sys.source(MX_PATHS$graphon_reference, envir = .GlobalEnv)

mx_gini <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x >= 0]
  if (!length(x) || sum(x) <= 0) return(0)
  x <- sort(x)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

mx_entropy_effective_n <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(0)
  p <- x / sum(x)
  exp(-sum(p * log(p)))
}

mx_weighted_reciprocity <- function(A) {
  A <- as_weight_matrix(A)
  diag(A) <- 0
  total <- sum(A)
  if (total <= 0) return(0)
  sum(pmin(A, t(A))) / total
}

mx_network_indicators <- function(x) {
  A <- as_weight_matrix(x)
  n <- nrow(A)
  diag_free <- A
  diag(diag_free) <- 0
  in_strength <- colSums(A)
  out_strength <- rowSums(A)
  total <- sum(A)
  top_n <- min(10L, n)
  top_share <- if (total > 0) {
    sum(sort(in_strength, decreasing = TRUE)[seq_len(top_n)]) / total
  } else {
    0
  }
  graph_binary <- igraph::graph_from_adjacency_matrix(
    (diag_free > 0) * 1,
    mode = "directed",
    weighted = NULL,
    diag = FALSE
  )
  transitivity <- suppressWarnings(
    igraph::transitivity(graph_binary, type = "global", isolates = "zero")
  )
  if (!is.finite(transitivity)) transitivity <- 0

  c(
    total_weight = total,
    edge_density = sum(diag_free > 0) / (n * (n - 1)),
    weighted_reciprocity = mx_weighted_reciprocity(A),
    in_strength_gini = mx_gini(in_strength),
    out_strength_gini = mx_gini(out_strength),
    effective_destinations = mx_entropy_effective_n(in_strength),
    top10_destination_share = top_share,
    binary_transitivity = transitivity
  )
}

mx_indicator_change_scores <- function(before, after) {
  before <- before[names(after)]
  if (anyNA(before) || anyNA(after)) stop("Indicator vectors are not aligned.")
  out <- abs(after - before)
  out["total_weight"] <- abs(
    log1p(after["total_weight"]) - log1p(before["total_weight"])
  )
  out["effective_destinations"] <- abs(
    log1p(after["effective_destinations"]) -
      log1p(before["effective_destinations"])
  )
  names(out) <- paste0("indicator_", names(out))
  out
}

mx_unit_mean_kernel <- function(M) {
  M <- as.matrix(M)
  M[!is.finite(M) | M < 0] <- 0
  scale_value <- mean(M)
  if (!is.finite(scale_value) || scale_value <= 0) return(matrix(0, nrow(M), ncol(M)))
  M / scale_value
}

mx_cut_distance_lower_bound <- function(M1, M2, samples = 1000L,
                                        restarts = 200L) {
  M1 <- as.matrix(M1)
  M2 <- as.matrix(M2)
  if (!identical(dim(M1), dim(M2))) stop("Kernels must have equal dimensions.")

  samples <- as.integer(samples)
  restarts <- as.integer(restarts)
  if (length(samples) != 1L || is.na(samples) || samples < 1L) {
    stop("samples must be one positive integer.", call. = FALSE)
  }
  if (length(restarts) != 1L || is.na(restarts) || restarts < 1L) {
    stop("restarts must be one positive integer.", call. = FALSE)
  }

  # The reference implementation uses alternating maximisation from
  # deterministic and random starts for both signs of the kernel difference.
  # Its previous random subset-pair estimate is retained as a floor. The
  # returned value is therefore a reproducible lower-bound estimate, not the
  # exact cut norm.
  cut_distance(
    M1,
    M2,
    num_samples = samples,
    cut_restarts = restarts
  )
}

mx_wasserstein_sorted <- function(M1, M2) {
  M1 <- as.numeric(M1)
  M2 <- as.numeric(M2)
  if (length(M1) != length(M2)) stop("Kernels must have equal sizes.")
  mean(abs(sort(M1) - sort(M2)))
}

mx_compact_graphon <- function(x) {
  list(
    W = x$W,
    fitted = x$fitted,
    K = x$K,
    n = x$n,
    node_labels = as.integer(x$node_labels),
    method = x$method,
    block_sizes = as.integer(tabulate(x$node_labels, nbins = x$K))
  )
}

mx_fit_graphon <- function(A, K, method = "block") {
  mx_compact_graphon(estimate_directed_graphon(A, K = K, method = method))
}

mx_pair_scores <- function(graphon_before, graphon_after, A_before, A_after,
                           cut_samples, cut_restarts) {
  F1 <- as.matrix(graphon_before$fitted)
  F2 <- as.matrix(graphon_after$fitted)
  if (!identical(dim(F1), dim(F2))) stop("Fitted kernels are not aligned.")
  S1 <- mx_unit_mean_kernel(F1)
  S2 <- mx_unit_mean_kernel(F2)
  A1 <- as_weight_matrix(A_before)
  A2 <- as_weight_matrix(A_after)
  AS1 <- mx_unit_mean_kernel(A1)
  AS2 <- mx_unit_mean_kernel(A2)
  n <- nrow(F1)

  c(
    graphon_frobenius_raw = norm(F1 - F2, type = "F") / n,
    graphon_cut_raw = mx_cut_distance_lower_bound(
      F1, F2, samples = cut_samples, restarts = cut_restarts
    ),
    graphon_wasserstein_raw = mx_wasserstein_sorted(F1, F2),
    graphon_frobenius_shape = norm(S1 - S2, type = "F") / n,
    graphon_cut_shape = mx_cut_distance_lower_bound(
      S1, S2, samples = cut_samples, restarts = cut_restarts
    ),
    graphon_wasserstein_shape = mx_wasserstein_sorted(S1, S2),
    singular_spectrum_raw = singular_spectrum_distance(A1, A2),
    singular_spectrum_shape = singular_spectrum_distance(AS1, AS2)
  )
}

mx_adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  observed <- sum(choose2(tab))
  row_total <- sum(choose2(rowSums(tab)))
  col_total <- sum(choose2(colSums(tab)))
  all_pairs <- choose2(sum(tab))
  if (all_pairs <= 0) return(NA_real_)
  expected <- row_total * col_total / all_pairs
  maximum <- 0.5 * (row_total + col_total)
  denominator <- maximum - expected
  if (abs(denominator) < .Machine$double.eps) return(1)
  (observed - expected) / denominator
}

mx_top_set <- function(x, n = 5L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(character())
  names(sort(x, decreasing = TRUE))[seq_len(min(n, length(x)))]
}

mx_jaccard <- function(a, b) {
  union_set <- union(a, b)
  if (!length(union_set)) return(NA_real_)
  length(intersect(a, b)) / length(union_set)
}

mx_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}

mx_dense_rank_desc <- function(x) {
  match(x, sort(unique(x), decreasing = TRUE))
}

mx_rank01 <- function(x) {
  if (length(x) <= 1L) return(rep(0.5, length(x)))
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

mx_role_profile <- function(graphon, A, node_metadata, year, core_definition) {
  labels <- as.integer(graphon$node_labels)
  K <- graphon$K
  sizes <- tabulate(labels, nbins = K)
  W <- graphon$W
  block_in <- as.numeric(t(W) %*% sizes) / sum(sizes)
  block_out <- as.numeric(W %*% sizes) / sum(sizes)
  destination_rank <- mx_dense_rank_desc(block_in)
  origin_rank <- mx_dense_rank_desc(block_out)
  core_idx <- match(node_metadata$nodeID, core_definition$nodeID)
  if (anyNA(core_idx)) stop("Destination-core alignment failed in ", year, ".")

  data.frame(
    year = as.integer(year),
    node_metadata,
    block = labels,
    block_size = sizes[labels],
    destination_role_rank = destination_rank[labels],
    origin_role_rank = origin_rank[labels],
    raw_in_strength = colSums(A),
    raw_out_strength = rowSums(A),
    fitted_in_mean = colMeans(graphon$fitted),
    fitted_out_mean = rowMeans(graphon$fitted),
    fitted_in_rank01 = mx_rank01(colMeans(graphon$fitted)),
    fitted_out_rank01 = mx_rank01(rowMeans(graphon$fitted)),
    destination_core_rank = core_definition$destination_core_rank[core_idx],
    is_destination_core = core_definition$is_destination_core[core_idx],
    is_regional_destination_core = core_definition$is_regional_destination_core[core_idx],
    stringsAsFactors = FALSE
  )
}

mx_same_group_rate <- function(labels) {
  labels <- labels[!is.na(labels)]
  if (length(labels) < 2L) return(NA_real_)
  same <- outer(labels, labels, FUN = "==")
  mean(same[upper.tri(same)])
}

mx_mean_pairwise_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  distance <- abs(outer(x, x, FUN = "-"))
  mean(distance[upper.tri(distance)])
}

mx_refinement_summary <- function(low_labels, high_labels, low_K, high_K, year) {
  tab <- table(low = low_labels, high = high_labels)
  high_sizes <- colSums(tab)
  dominant <- apply(tab, 2L, max)
  shares <- dominant / high_sizes
  data.frame(
    year = as.integer(year),
    low_K = as.integer(low_K),
    high_K = as.integer(high_K),
    weighted_refinement_consistency = sum(dominant) / sum(tab),
    minimum_high_block_parent_share = min(shares),
    median_high_block_parent_share = stats::median(shares),
    stringsAsFactors = FALSE
  )
}

mx_log("Analytical helpers loaded.")



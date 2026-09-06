# =============================================================================
# DIRECTED WEIGHTED GRAPHON / BLOCK-KERNEL ANALYSIS
# Revision-safe implementation for the Mob project
# =============================================================================
#
# Design principles for the revision:
#   1. The manuscript's primary temporal analysis uses fixed K=13. The
#      snapshot-specific adaptive K_t series is diagnostic only. NOTE
#      (2026-08): the block-number heuristic previously returned the
#      floor(sqrt(n)) ceiling whenever a yearly matrix was rank deficient,
#      because a numerically zero singular value dominated the relative-drop
#      ratios. estimate_number_of_blocks() now removes numerical zeros before
#      applying the unchanged relative-drop rule.
#   2. The original weighted scale is retained in the main graphon estimates.
#   3. Distances between graphon estimates are computed on node-aligned fitted
#      n x n kernel matrices. Therefore arbitrary block-label permutations and
#      differing K_t values do not require zero-padding of K x K matrices.
#   4. Directed spectral comparison is based on singular values of the original
#      weighted adjacency matrices and is reported separately as an adjacency
#      benchmark, not as a graphon-estimator distance.
#   5. Fixed-K estimation supports the primary and sensitivity analyses.
#   6. Missing packages cause an explicit error; this file never installs them.
# =============================================================================

required_packages <- c("igraph", "Matrix", "RSpectra", "ggplot2", "reshape2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

as_weight_matrix <- function(x) {
  if (inherits(x, "igraph")) {
    A <- as.matrix(igraph::as_adjacency_matrix(x, attr = "weight", sparse = FALSE))
  } else {
    A <- as.matrix(x)
  }
  if (nrow(A) != ncol(A)) stop("Adjacency matrix must be square.")
  storage.mode(A) <- "double"
  A[!is.finite(A)] <- 0
  A[A < 0] <- 0
  A
}

safe_scale <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

block_features <- function(A) {
  A <- as_weight_matrix(A)
  cbind(
    out_strength = safe_scale(rowSums(A)),
    in_strength  = safe_scale(colSums(A)),
    out_degree   = safe_scale(rowSums(A > 0)),
    in_degree    = safe_scale(colSums(A > 0))
  )
}

n_distinct_rows <- function(M, digits = 12L) {
  M <- as.matrix(M)
  if (!nrow(M)) return(0L)
  M[!is.finite(M)] <- 0
  as.integer(nrow(unique(round(M, digits = digits))))
}

estimate_number_of_blocks <- function(A) {
  A <- as_weight_matrix(A)
  n <- nrow(A)
  if (n <= 1L) return(1L)

  # Preserve the original manuscript's elbow search range: inspect up to the
  # first 20 singular values, and apply sqrt(n) only to the selected K.
  max_k <- min(20L, n - 1L)
  if (max_k < 2L) return(1L)

  # FIX (2026-08): use the dense, deterministic svd() here. The iterative
  # RSpectra solver has a noise floor around 1e-9 relative, which is large
  # enough to hide a numerically zero singular value and therefore defeats the
  # rank check below. For n <= 200 the dense decomposition costs milliseconds.
  d <- tryCatch({
    if (n > 200L) {
      RSpectra::svds(A, k = max_k)$d
    } else {
      svd(A, nu = 0, nv = 0)$d[seq_len(max_k)]
    }
  }, error = function(e) {
    svd(A, nu = 0, nv = 0)$d[seq_len(max_k)]
  })

  d <- as.numeric(d)
  d <- d[is.finite(d) & d >= 0]
  if (length(d) < 2L || max(d) <= .Machine$double.eps) return(1L)
  if (length(d) < 3L) return(2L)

  # FIX (2026-08): drop numerically zero singular values before forming the
  # relative-drop ratios. Several yearly matrices have rank 19, so d[20] is
  # ~1e-10. With the previous formulation d[19] / (d[20] + 1e-10) exploded to
  # ~1e11 and which.max() always selected the last index, forcing K to the
  # floor(sqrt(n)) ceiling. K then reflected rank deficiency rather than
  # structural differentiation. The relative-drop rule itself is unchanged.
  # Relative tolerance. On these data the separation is unambiguous: in
  # rank-deficient years the last singular value is ~1e-15 of the largest,
  # while in full-rank years it is ~1e-4. A 1e-8 relative cut sits safely
  # between the two.
  tol <- max(d) * 1e-8
  d <- d[d > tol]
  if (length(d) < 3L) return(2L)

  ratios <- d[-length(d)] / d[-1L]
  candidate <- as.integer(which.max(ratios) + 1L)
  candidate <- max(2L, min(candidate, floor(sqrt(n))))

  # Technical estimability guard only: an automatically selected K cannot
  # exceed the number of distinct clustering feature profiles. This does not
  # alter AGGR snapshots when the original candidate is estimable.
  feasible <- n_distinct_rows(block_features(A))
  if (feasible <= 1L) return(1L)
  as.integer(min(candidate, feasible))
}

run_kmeans_checked <- function(features, K, context = "graphon estimator") {
  features <- as.matrix(features)
  features[!is.finite(features)] <- 0
  distinct_n <- n_distinct_rows(features)

  if (K < 1L) stop("K must be at least 1.")
  if (K > distinct_n) {
    stop(
      context, ": requested K=", K,
      " but only ", distinct_n, " distinct feature vectors are available."
    )
  }
  if (K == 1L) return(rep.int(1L, nrow(features)))

  set.seed(42)
  stats::kmeans(features, centers = K, nstart = 25, iter.max = 100)$cluster
}

block_matrix_from_labels <- function(A, labels, K) {
  W <- matrix(0, K, K)
  for (i in seq_len(K)) {
    ii <- which(labels == i)
    for (j in seq_len(K)) {
      jj <- which(labels == j)
      if (length(ii) > 0L && length(jj) > 0L) {
        W[i, j] <- mean(A[ii, jj, drop = FALSE])
      }
    }
  }
  W
}

estimate_graphon_block <- function(A, K) {
  features <- block_features(A)
  labels <- run_kmeans_checked(features, K, context = "block estimator")
  list(W = block_matrix_from_labels(A, labels, K), labels = labels)
}

estimate_graphon_smooth <- function(A, K) {
  n <- nrow(A)
  if (K == 1L) {
    labels <- rep.int(1L, n)
    return(list(W = matrix(mean(A), 1L, 1L), labels = labels))
  }

  k_use <- min(max(K, 2L), n - 1L)
  sv <- if (n > 100L) RSpectra::svds(A, k = k_use) else svd(A)
  U <- sv$u[, seq_len(k_use), drop = FALSE]
  V <- sv$v[, seq_len(k_use), drop = FALSE]
  d <- sv$d[seq_len(k_use)]

  A_smooth <- U %*% diag(d, nrow = k_use) %*% t(V)
  A_smooth[!is.finite(A_smooth) | A_smooth < 0] <- 0

  features <- cbind(U, V)
  labels <- run_kmeans_checked(features, K, context = "smooth estimator")
  list(W = block_matrix_from_labels(A_smooth, labels, K), labels = labels)
}

estimate_graphon_spectral <- function(A, K) {
  n <- nrow(A)
  if (K == 1L) {
    labels <- rep.int(1L, n)
    return(list(W = matrix(mean(A), 1L, 1L), labels = labels))
  }

  k_use <- min(max(K, 2L), n - 1L)
  sv <- if (n > 100L) RSpectra::svds(A, k = k_use) else svd(A)
  U <- sv$u[, seq_len(k_use), drop = FALSE]
  V <- sv$v[, seq_len(k_use), drop = FALSE]
  d <- sv$d[seq_len(k_use)]

  embedding <- cbind(
    U %*% diag(sqrt(pmax(d, 0)), nrow = k_use),
    V %*% diag(sqrt(pmax(d, 0)), nrow = k_use)
  )
  embedding[!is.finite(embedding)] <- 0
  labels <- run_kmeans_checked(embedding, K, context = "spectral estimator")
  list(W = block_matrix_from_labels(A, labels, K), labels = labels)
}

fitted_kernel_matrix <- function(W, labels) {
  W[labels, labels, drop = FALSE]
}

estimate_directed_graphon <- function(A, K = NULL,
                                      method = c("block", "smooth", "spectral")) {
  method <- match.arg(method)
  A <- as_weight_matrix(A)
  candidate_K <- if (is.null(K)) estimate_number_of_blocks(A) else as.integer(K)

  if (length(candidate_K) != 1L || !is.finite(candidate_K)) {
    stop("K must be NULL or one finite integer.")
  }
  max_K <- max(1L, floor(sqrt(nrow(A))))
  if (candidate_K < 1L || candidate_K > max_K) {
    stop("K must be between 1 and floor(sqrt(n)).")
  }

  result <- switch(
    method,
    block    = estimate_graphon_block(A, candidate_K),
    smooth   = estimate_graphon_smooth(A, candidate_K),
    spectral = estimate_graphon_spectral(A, candidate_K)
  )

  out <- list(
    W = result$W,
    fitted = fitted_kernel_matrix(result$W, result$labels),
    K = candidate_K,
    n = nrow(A),
    node_labels = result$labels,
    method = method,
    block_sizes = table(result$labels),
    original_matrix = A,
    raw_total_weight = sum(A),
    normalization = "none",
    K_selection = if (is.null(K)) "adaptive_snapshot" else "fixed_or_supplied"
  )
  class(out) <- "directed_graphon"
  out
}

get_block_membership <- function(graphon, node_names = NULL,
                                 output_format = c("list", "data.frame")) {
  if (!inherits(graphon, "directed_graphon")) {
    stop("Input must be a directed_graphon object.")
  }
  output_format <- match.arg(output_format)
  ids <- if (is.null(node_names)) seq_len(graphon$n) else node_names
  if (length(ids) != graphon$n) stop("node_names has incorrect length.")

  if (output_format == "data.frame") {
    out <- data.frame(node = ids, block = graphon$node_labels, stringsAsFactors = FALSE)
    out <- out[order(out$block, out$node), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  out <- split(ids, factor(graphon$node_labels, levels = seq_len(graphon$K)))
  names(out) <- paste0("Block_", seq_len(graphon$K))
  out
}

top_singular_values <- function(A, k) {
  A <- as_weight_matrix(A)
  if (k < 1L) return(numeric())
  tryCatch(
    {
      if (nrow(A) > 50L) {
        as.numeric(RSpectra::svds(A, k = k)$d)
      } else {
        as.numeric(svd(A, nu = 0, nv = 0)$d[seq_len(k)])
      }
    },
    error = function(e) as.numeric(svd(A, nu = 0, nv = 0)$d[seq_len(k)])
  )
}

singular_spectrum_distance <- function(A1, A2, k = 10L) {
  A1 <- as_weight_matrix(A1)
  A2 <- as_weight_matrix(A2)
  if (!identical(dim(A1), dim(A2))) {
    stop("Adjacency matrices must have equal dimensions.")
  }

  k <- min(as.integer(k), nrow(A1) - 1L, nrow(A2) - 1L)
  if (k < 1L) return(0)

  s1 <- top_singular_values(A1, k)
  s2 <- top_singular_values(A2, k)
  sqrt(sum((s1 - s2)^2)) / sqrt(k)
}

spectral_distance <- function(A1, A2, k = 10L) {
  singular_spectrum_distance(A1, A2, k = k)
}

as_fitted_matrix <- function(x) {
  if (inherits(x, "directed_graphon")) return(x$fitted)
  as.matrix(x)
}

frobenius_distance <- function(M1, M2, A1 = NULL, A2 = NULL) {
  if (is.null(M1) || is.null(M2)) {
    if (is.null(A1) || is.null(A2)) return(NA_real_)
    M1 <- as_weight_matrix(A1)
    M2 <- as_weight_matrix(A2)
  } else {
    M1 <- as_fitted_matrix(M1)
    M2 <- as_fitted_matrix(M2)
  }
  if (!identical(dim(M1), dim(M2))) stop("Fitted kernels must have equal dimensions.")
  norm(M1 - M2, type = "F") / nrow(M1)
}

cut_distance <- function(M1, M2, A1 = NULL, A2 = NULL, num_samples = 1000L,
                         cut_restarts = 200L) {
  if (is.null(M1) || is.null(M2)) {
    if (is.null(A1) || is.null(A2)) return(NA_real_)
    M1 <- as_weight_matrix(A1)
    M2 <- as_weight_matrix(A2)
  } else {
    M1 <- as_fitted_matrix(M1)
    M2 <- as_fitted_matrix(M2)
  }
  if (!identical(dim(M1), dim(M2))) stop("Fitted kernels must have equal dimensions.")

  D <- M1 - M2
  n <- nrow(D)
  set.seed(42)

  # FIX (2026-08): the cut norm is a maximum over all subset pairs (S, T).
  # Drawing num_samples random pairs from a 2^n x 2^n space is a very loose
  # lower bound: on these data it disagreed with a proper search at rank
  # correlation 0.34 and changed which year pair ranked first.
  #
  # Alternating maximisation is used instead. For fixed S the optimal T is
  # available in closed form (columns whose column sum over S is positive)
  # and vice versa, so the iteration increases the objective monotonically
  # and terminates. Deterministic and random restarts are run for both +D
  # and -D. The original random-sampling result is retained as a floor, so
  # the returned value is never below the previous implementation.
  #
  # NOTE: this is a tighter lower bound, not an exact value. Exact
  # computation of the cut norm is NP-hard. Report it as
  # "alternating maximisation with N random restarts".
  ascend <- function(M, S) {
    for (it in seq_len(100L)) {
      cs <- colSums(M[S, , drop = FALSE])
      Tn <- cs > 0
      if (!any(Tn)) Tn <- cs >= max(cs)
      rs <- rowSums(M[, Tn, drop = FALSE])
      Sn <- rs > 0
      if (!any(Sn)) Sn <- rs >= max(rs)
      if (identical(Sn, S)) return(sum(M[Sn, Tn, drop = FALSE]))
      S <- Sn
    }
    cs <- colSums(M[S, , drop = FALSE])
    Tn <- cs > 0
    if (!any(Tn)) Tn <- cs >= max(cs)
    sum(M[S, Tn, drop = FALSE])
  }

  max_cut <- 0
  starts_det <- list(rep(TRUE, n), rowSums(D) > 0, rowSums(D) < 0)
  for (M in list(D, -D)) {
    for (S0 in starts_det) {
      if (any(S0)) max_cut <- max(max_cut, ascend(M, S0))
    }
    for (r in seq_len(as.integer(cut_restarts))) {
      S0 <- sample(c(TRUE, FALSE), n, replace = TRUE)
      if (any(S0)) max_cut <- max(max_cut, ascend(M, S0))
    }
  }

  # Floor from the previous random-sampling estimator.
  for (i in seq_len(as.integer(num_samples))) {
    S <- sample(c(TRUE, FALSE), n, replace = TRUE)
    T <- sample(c(TRUE, FALSE), n, replace = TRUE)
    max_cut <- max(max_cut, abs(sum(D[S, T, drop = FALSE])))
  }
  max_cut / (n * n)
}

wasserstein_distance <- function(M1, M2) {
  M1 <- as_fitted_matrix(M1)
  M2 <- as_fitted_matrix(M2)
  if (length(M1) != length(M2)) stop("Fitted kernels must have equal sizes.")
  mean(abs(sort(as.numeric(M1)) - sort(as.numeric(M2))))
}

js_divergence <- function(M1, M2) {
  M1 <- as_fitted_matrix(M1)
  M2 <- as_fitted_matrix(M2)
  p <- pmax(as.numeric(M1), 0)
  q <- pmax(as.numeric(M2), 0)
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)

  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  kl <- function(a, b) {
    keep <- a > 0 & b > 0
    sum(a[keep] * log(a[keep] / b[keep]))
  }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

tv_distance <- function(M1, M2) {
  M1 <- as_fitted_matrix(M1)
  M2 <- as_fitted_matrix(M2)
  p <- pmax(as.numeric(M1), 0)
  q <- pmax(as.numeric(M2), 0)
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)

  p <- p / sum(p)
  q <- q / sum(q)
  0.5 * sum(abs(p - q))
}

compute_graphon_distance <- function(graphon1, graphon2,
                                     method = c("all", "frobenius", "cut", "wasserstein",
                                                "js", "tv", "singular_spectrum", "spectral")) {
  method <- match.arg(method)
  F1 <- graphon1$fitted
  F2 <- graphon2$fitted
  A1 <- graphon1$original_matrix
  A2 <- graphon2$original_matrix

  if (!identical(dim(F1), dim(F2))) {
    stop("Graphons must refer to the same node set for node-aligned comparison.")
  }

  if (method == "all") {
    return(list(
      singular_spectrum = singular_spectrum_distance(A1, A2),
      frobenius = frobenius_distance(F1, F2),
      cut = cut_distance(F1, F2),
      wasserstein = wasserstein_distance(F1, F2),
      jensen_shannon = js_divergence(F1, F2),
      total_variation = tv_distance(F1, F2)
    ))
  }
  if (method %in% c("singular_spectrum", "spectral")) return(singular_spectrum_distance(A1, A2))
  if (method == "frobenius") return(frobenius_distance(F1, F2))
  if (method == "cut") return(cut_distance(F1, F2))
  if (method == "wasserstein") return(wasserstein_distance(F1, F2))
  if (method == "js") return(js_divergence(F1, F2))
  if (method == "tv") return(tv_distance(F1, F2))
}

graphon_distance <- compute_graphon_distance

resolve_K_schedule <- function(network_list, K = NULL) {
  n_net <- length(network_list)
  if (is.null(K)) {
    out <- vapply(network_list, function(x) estimate_number_of_blocks(as_weight_matrix(x)), integer(1))
    return(out)
  }

  K <- as.integer(K)
  if (length(K) == 1L) return(rep.int(K, n_net))
  if (length(K) != n_net) {
    stop("K must be NULL, one integer, or one integer per network.")
  }
  K
}

compare_graphons <- function(network_list, K = NULL,
                             method = c("block", "smooth", "spectral")) {
  method <- match.arg(method)
  if (length(network_list) < 2L) stop("At least two networks are required.")

  network_names <- names(network_list)
  if (is.null(network_names)) network_names <- as.character(seq_along(network_list))

  dims <- vapply(network_list, function(x) nrow(as_weight_matrix(x)), integer(1))
  if (length(unique(dims)) != 1L) {
    stop("All networks in a comparison set must contain the same number of nodes.")
  }

  K_schedule <- resolve_K_schedule(network_list, K = K)
  names(K_schedule) <- network_names

  graphons <- lapply(seq_along(network_list), function(i) {
    estimate_directed_graphon(network_list[[i]], K = K_schedule[i], method = method)
  })
  names(graphons) <- network_names

  metric_names <- c(
    "singular_spectrum", "frobenius", "cut", "wasserstein",
    "jensen_shannon", "total_variation"
  )
  distances <- setNames(lapply(metric_names, function(x) {
    matrix(0, length(graphons), length(graphons),
           dimnames = list(network_names, network_names))
  }), metric_names)

  for (i in seq_len(length(graphons) - 1L)) {
    for (j in seq.int(i + 1L, length(graphons))) {
      d <- compute_graphon_distance(graphons[[i]], graphons[[j]], method = "all")
      for (nm in metric_names) {
        distances[[nm]][i, j] <- distances[[nm]][j, i] <- d[[nm]]
      }
    }
  }

  distances$spectral <- distances$singular_spectrum

  list(
    graphons = graphons,
    distances = distances,
    network_names = network_names,
    K = K_schedule,
    K_by_network = K_schedule,
    method = method,
    normalization = "none",
    distance_domain = "node_aligned_fitted_kernel",
    benchmark = "singular_spectrum_on_original_weighted_adjacency"
  )
}

plot_graphon <- function(graphon, title = NULL) {
  if (!inherits(graphon, "directed_graphon")) {
    stop("graphon must be a directed_graphon object.")
  }
  df <- reshape2::melt(graphon$W)
  names(df) <- c("from_block", "to_block", "weight")
  ggplot2::ggplot(df, ggplot2::aes(x = to_block, y = from_block, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = "Destination block", y = "Origin block", fill = "Mean weight",
      title = if (is.null(title)) paste0("Directed weighted block kernel (K=", graphon$K, ")") else title
    ) +
    ggplot2::theme_minimal()
}

plot_distance_heatmap <- function(result, metric = "frobenius") {
  D <- result$distances[[metric]]
  if (is.null(D)) stop("Unknown distance metric: ", metric)
  df <- reshape2::melt(D)
  names(df) <- c("network_1", "network_2", "distance")
  ggplot2::ggplot(df, ggplot2::aes(x = network_2, y = network_1, fill = distance)) +
    ggplot2::geom_tile() +
    ggplot2::labs(x = NULL, y = NULL, fill = metric) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}

generate_pastel_palette <- function(n) {
  if (n <= 0L) return(character())
  hues <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(h = hues, c = 45, l = 75)
}

mix_colors_additive <- function(colors) {
  if (length(colors) == 0L) return("gray70")
  rgb_mtx <- grDevices::col2rgb(colors) / 255
  m <- pmin(1, rowMeans(rgb_mtx))
  grDevices::rgb(m[1], m[2], m[3])
}

darken_color <- function(color, factor = 0.75) {
  rgb <- grDevices::col2rgb(color) * factor
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
}

annotate_graph_with_graphon <- function(g,
                                        K = NULL,
                                        method = "block",
                                        block_colors = NULL,
                                        default_color = "gray70",
                                        log_transform = c("none", "log", "loglog"),
                                        min_width = 0.5,
                                        max_width = 6) {
  log_transform <- match.arg(log_transform)
  if (!inherits(g, "igraph")) stop("g must be an igraph object.")

  A <- as_weight_matrix(g)
  graphon_obj <- estimate_directed_graphon(A, K = K, method = method)
  member <- graphon_obj$node_labels

  g <- igraph::set_vertex_attr(g, "graphonmember", value = member)
  out_strength <- rowSums(A)
  in_strength <- colSums(A)
  g <- igraph::set_vertex_attr(
    g, "origmember", value = ifelse(out_strength > 0, member, NA_integer_)
  )
  g <- igraph::set_vertex_attr(
    g, "destmember", value = ifelse(in_strength > 0, member, NA_integer_)
  )

  if (is.null(block_colors)) block_colors <- generate_pastel_palette(graphon_obj$K)
  if (length(block_colors) < graphon_obj$K) {
    block_colors <- rep(block_colors, length.out = graphon_obj$K)
  }
  if (length(block_colors) > 0L) {
    g <- igraph::set_vertex_attr(g, "color", value = block_colors[member])
  } else {
    g <- igraph::set_vertex_attr(g, "color", value = rep(default_color, igraph::vcount(g)))
  }

  if (igraph::ecount(g) > 0L) {
    ends <- igraph::ends(g, igraph::E(g), names = FALSE)
    edge_cols <- vapply(seq_len(nrow(ends)), function(i) {
      mix_colors_additive(block_colors[member[ends[i, ]]])
    }, character(1))
    g <- igraph::set_edge_attr(g, "color", value = edge_cols)

    w <- if ("weight" %in% igraph::edge_attr_names(g)) {
      as.numeric(igraph::edge_attr(g, "weight"))
    } else {
      rep(1, igraph::ecount(g))
    }
    if (log_transform == "log") w <- log1p(pmax(w, 0))
    if (log_transform == "loglog") w <- log1p(log1p(pmax(w, 0)))

    finite_w <- w[is.finite(w)]
    if (length(finite_w) == 0L || length(unique(finite_w)) <= 1L) {
      widths <- rep((min_width + max_width) / 2, length(w))
    } else {
      r <- range(finite_w)
      widths <- min_width + (w - r[1]) / (r[2] - r[1]) * (max_width - min_width)
      widths[!is.finite(widths)] <- min_width
    }
    g <- igraph::set_edge_attr(g, "width", value = widths)
  }

  attr(g, "graphon") <- graphon_obj
  g
}




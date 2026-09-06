# =============================================================================
# run_all_standalone.R
# Standalone reproducibility pipeline for the Mob project
# =============================================================================
#
# Usage from the project directory:
#   source("run_all_standalone.R")
# or:
#   Rscript run_all_standalone.R
#
# Required project inputs:
#   data/
#   run_all_standalone.R
#
# For manuscript rendering after the run:
#   Mob_final.Rmd
#   bibliography-final.bib
#
# The script does not install packages and does not source project-local R code.
# The non-CRAN tsnda package is NOT required: the author's vgprops() and
# percolate() implementations are embedded below. The pipeline creates the
# calculation objects, reviewer outputs, figure-data tables and vector-PDF
# figures consumed by Mob_final.Rmd.
#
# Reproducibility controls:
#   MOB_NLAC_REPS              default 10000 (legacy definition)
#   MOB_SKIP_LEGACY_GRAPHONS   default 0; set 1 only for diagnostics
#   MOB_REFRESH_RENDER_HELPERS default 1; set 0 only to retain existing helpers
# =============================================================================

options(stringsAsFactors = FALSE)

.locate_mob_root <- function() {
  candidates <- character()
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    candidates <- c(candidates, dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    candidates <- c(candidates, dirname(normalizePath(sub("^--file=", "", file_arg[1L]),
                                                       winslash = "/", mustWork = FALSE)))
  }
  candidates <- unique(c(candidates, getwd()))
  ok <- vapply(candidates, function(x) file.exists(file.path(x, "data", "nodes.txt")), logical(1))
  if (!any(ok)) {
    stop("Cannot locate the Mob project root. Run/source this file from a directory containing data/nodes.txt.")
  }
  normalizePath(candidates[which(ok)[1L]], winslash = "/", mustWork = TRUE)
}

ROOT <- .locate_mob_root()
setwd(ROOT)

SEED <- 20260811L
set.seed(SEED)
NLAC_REPS <- suppressWarnings(as.integer(Sys.getenv("MOB_NLAC_REPS", "10000")))
if (!is.finite(NLAC_REPS) || NLAC_REPS < 1L) stop("MOB_NLAC_REPS must be a positive integer.")
SKIP_LEGACY_GRAPHONS <- identical(Sys.getenv("MOB_SKIP_LEGACY_GRAPHONS", "0"), "1")
REFRESH_RENDER_HELPERS <- identical(Sys.getenv("MOB_REFRESH_RENDER_HELPERS", "1"), "1")

dir.create("calcs", recursive = TRUE, showWarnings = FALSE)
dir.create("calcs/revision", recursive = TRUE, showWarnings = FALSE)
dir.create("output", recursive = TRUE, showWarnings = FALSE)
dir.create("output/revision", recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- "calcs/run_all_standalone_runtime.log"
writeLines(character(), LOG_FILE)

log_step <- function(...) {
  msg <- paste0(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
  invisible(msg)
}

required_packages <- c(
  "igraph", "brainGraph", "Matrix", "RSpectra", "irlba",
  "ggplot2", "reshape2", "dplyr", "tidyr", "openxlsx", "nda"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before running; this script never installs packages automatically."
  )
}

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

assert_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required file(s): ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
assert_files("data/nodes.txt")

safe_scalar <- function(expr, context) {
  z <- tryCatch(force(expr), error = function(e) {
    stop(context, ": ", conditionMessage(e), call. = FALSE)
  })
  z <- as.numeric(z)
  if (!length(z)) stop(context, ": empty result.", call. = FALSE)
  z[1L]
}

write_embedded_helper <- function(path, text) {
  if (!file.exists(path) || REFRESH_RENDER_HELPERS) {
    writeLines(text, path, useBytes = TRUE)
    log_step("Wrote render-compatibility helper: ", path)
  }
  invisible(path)
}

GRAPHON_HELPER_SOURCE <- "# =============================================================================\n# DIRECTED WEIGHTED GRAPHON / BLOCK-KERNEL ANALYSIS\n# Revision-safe implementation for the Mob project\n# =============================================================================\n#\n# Design principles for the revision:\n#   1. The manuscript's primary temporal analysis uses fixed K=13. The\n#      snapshot-specific adaptive K_t series is diagnostic only. NOTE\n#      (2026-08): the block-number heuristic previously returned the\n#      floor(sqrt(n)) ceiling whenever a yearly matrix was rank deficient,\n#      because a numerically zero singular value dominated the relative-drop\n#      ratios. estimate_number_of_blocks() now removes numerical zeros before\n#      applying the unchanged relative-drop rule.\n#   2. The original weighted scale is retained in the main graphon estimates.\n#   3. Distances between graphon estimates are computed on node-aligned fitted\n#      n x n kernel matrices. Therefore arbitrary block-label permutations and\n#      differing K_t values do not require zero-padding of K x K matrices.\n#   4. Directed spectral comparison is based on singular values of the original\n#      weighted adjacency matrices and is reported separately as an adjacency\n#      benchmark, not as a graphon-estimator distance.\n#   5. Fixed-K estimation supports the primary and sensitivity analyses.\n#   6. Missing packages cause an explicit error; this file never installs them.\n# =============================================================================\n\nrequired_packages <- c(\"igraph\", \"Matrix\", \"RSpectra\", \"ggplot2\", \"reshape2\")\nmissing_packages <- required_packages[\n  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)\n]\nif (length(missing_packages) > 0L) {\n  stop(\"Missing required R package(s): \", paste(missing_packages, collapse = \", \"))\n}\n\nas_weight_matrix <- function(x) {\n  if (inherits(x, \"igraph\")) {\n    A <- as.matrix(igraph::as_adjacency_matrix(x, attr = \"weight\", sparse = FALSE))\n  } else {\n    A <- as.matrix(x)\n  }\n  if (nrow(A) != ncol(A)) stop(\"Adjacency matrix must be square.\")\n  storage.mode(A) <- \"double\"\n  A[!is.finite(A)] <- 0\n  A[A < 0] <- 0\n  A\n}\n\nsafe_scale <- function(x) {\n  s <- stats::sd(x)\n  if (!is.finite(s) || s == 0) return(rep(0, length(x)))\n  as.numeric(scale(x))\n}\n\nblock_features <- function(A) {\n  A <- as_weight_matrix(A)\n  cbind(\n    out_strength = safe_scale(rowSums(A)),\n    in_strength  = safe_scale(colSums(A)),\n    out_degree   = safe_scale(rowSums(A > 0)),\n    in_degree    = safe_scale(colSums(A > 0))\n  )\n}\n\nn_distinct_rows <- function(M, digits = 12L) {\n  M <- as.matrix(M)\n  if (!nrow(M)) return(0L)\n  M[!is.finite(M)] <- 0\n  as.integer(nrow(unique(round(M, digits = digits))))\n}\n\nestimate_number_of_blocks <- function(A) {\n  A <- as_weight_matrix(A)\n  n <- nrow(A)\n  if (n <= 1L) return(1L)\n\n  # Preserve the original manuscript's elbow search range: inspect up to the\n  # first 20 singular values, and apply sqrt(n) only to the selected K.\n  max_k <- min(20L, n - 1L)\n  if (max_k < 2L) return(1L)\n\n  # FIX (2026-08): use the dense, deterministic svd() here. The iterative\n  # RSpectra solver has a noise floor around 1e-9 relative, which is large\n  # enough to hide a numerically zero singular value and therefore defeats the\n  # rank check below. For n <= 200 the dense decomposition costs milliseconds.\n  d <- tryCatch({\n    if (n > 200L) {\n      RSpectra::svds(A, k = max_k)$d\n    } else {\n      svd(A, nu = 0, nv = 0)$d[seq_len(max_k)]\n    }\n  }, error = function(e) {\n    svd(A, nu = 0, nv = 0)$d[seq_len(max_k)]\n  })\n\n  d <- as.numeric(d)\n  d <- d[is.finite(d) & d >= 0]\n  if (length(d) < 2L || max(d) <= .Machine$double.eps) return(1L)\n  if (length(d) < 3L) return(2L)\n\n  # FIX (2026-08): drop numerically zero singular values before forming the\n  # relative-drop ratios. Several yearly matrices have rank 19, so d[20] is\n  # ~1e-10. With the previous formulation d[19] / (d[20] + 1e-10) exploded to\n  # ~1e11 and which.max() always selected the last index, forcing K to the\n  # floor(sqrt(n)) ceiling. K then reflected rank deficiency rather than\n  # structural differentiation. The relative-drop rule itself is unchanged.\n  # Relative tolerance. On these data the separation is unambiguous: in\n  # rank-deficient years the last singular value is ~1e-15 of the largest,\n  # while in full-rank years it is ~1e-4. A 1e-8 relative cut sits safely\n  # between the two.\n  tol <- max(d) * 1e-8\n  d <- d[d > tol]\n  if (length(d) < 3L) return(2L)\n\n  ratios <- d[-length(d)] / d[-1L]\n  candidate <- as.integer(which.max(ratios) + 1L)\n  candidate <- max(2L, min(candidate, floor(sqrt(n))))\n\n  # Technical estimability guard only: an automatically selected K cannot\n  # exceed the number of distinct clustering feature profiles. This does not\n  # alter AGGR snapshots when the original candidate is estimable.\n  feasible <- n_distinct_rows(block_features(A))\n  if (feasible <= 1L) return(1L)\n  as.integer(min(candidate, feasible))\n}\n\nrun_kmeans_checked <- function(features, K, context = \"graphon estimator\") {\n  features <- as.matrix(features)\n  features[!is.finite(features)] <- 0\n  distinct_n <- n_distinct_rows(features)\n\n  if (K < 1L) stop(\"K must be at least 1.\")\n  if (K > distinct_n) {\n    stop(\n      context, \": requested K=\", K,\n      \" but only \", distinct_n, \" distinct feature vectors are available.\"\n    )\n  }\n  if (K == 1L) return(rep.int(1L, nrow(features)))\n\n  set.seed(42)\n  stats::kmeans(features, centers = K, nstart = 25, iter.max = 100)$cluster\n}\n\nblock_matrix_from_labels <- function(A, labels, K) {\n  W <- matrix(0, K, K)\n  for (i in seq_len(K)) {\n    ii <- which(labels == i)\n    for (j in seq_len(K)) {\n      jj <- which(labels == j)\n      if (length(ii) > 0L && length(jj) > 0L) {\n        W[i, j] <- mean(A[ii, jj, drop = FALSE])\n      }\n    }\n  }\n  W\n}\n\nestimate_graphon_block <- function(A, K) {\n  features <- block_features(A)\n  labels <- run_kmeans_checked(features, K, context = \"block estimator\")\n  list(W = block_matrix_from_labels(A, labels, K), labels = labels)\n}\n\nestimate_graphon_smooth <- function(A, K) {\n  n <- nrow(A)\n  if (K == 1L) {\n    labels <- rep.int(1L, n)\n    return(list(W = matrix(mean(A), 1L, 1L), labels = labels))\n  }\n\n  k_use <- min(max(K, 2L), n - 1L)\n  sv <- if (n > 100L) RSpectra::svds(A, k = k_use) else svd(A)\n  U <- sv$u[, seq_len(k_use), drop = FALSE]\n  V <- sv$v[, seq_len(k_use), drop = FALSE]\n  d <- sv$d[seq_len(k_use)]\n\n  A_smooth <- U %*% diag(d, nrow = k_use) %*% t(V)\n  A_smooth[!is.finite(A_smooth) | A_smooth < 0] <- 0\n\n  features <- cbind(U, V)\n  labels <- run_kmeans_checked(features, K, context = \"smooth estimator\")\n  list(W = block_matrix_from_labels(A_smooth, labels, K), labels = labels)\n}\n\nestimate_graphon_spectral <- function(A, K) {\n  n <- nrow(A)\n  if (K == 1L) {\n    labels <- rep.int(1L, n)\n    return(list(W = matrix(mean(A), 1L, 1L), labels = labels))\n  }\n\n  k_use <- min(max(K, 2L), n - 1L)\n  sv <- if (n > 100L) RSpectra::svds(A, k = k_use) else svd(A)\n  U <- sv$u[, seq_len(k_use), drop = FALSE]\n  V <- sv$v[, seq_len(k_use), drop = FALSE]\n  d <- sv$d[seq_len(k_use)]\n\n  embedding <- cbind(\n    U %*% diag(sqrt(pmax(d, 0)), nrow = k_use),\n    V %*% diag(sqrt(pmax(d, 0)), nrow = k_use)\n  )\n  embedding[!is.finite(embedding)] <- 0\n  labels <- run_kmeans_checked(embedding, K, context = \"spectral estimator\")\n  list(W = block_matrix_from_labels(A, labels, K), labels = labels)\n}\n\nfitted_kernel_matrix <- function(W, labels) {\n  W[labels, labels, drop = FALSE]\n}\n\nestimate_directed_graphon <- function(A, K = NULL,\n                                      method = c(\"block\", \"smooth\", \"spectral\")) {\n  method <- match.arg(method)\n  A <- as_weight_matrix(A)\n  candidate_K <- if (is.null(K)) estimate_number_of_blocks(A) else as.integer(K)\n\n  if (length(candidate_K) != 1L || !is.finite(candidate_K)) {\n    stop(\"K must be NULL or one finite integer.\")\n  }\n  max_K <- max(1L, floor(sqrt(nrow(A))))\n  if (candidate_K < 1L || candidate_K > max_K) {\n    stop(\"K must be between 1 and floor(sqrt(n)).\")\n  }\n\n  result <- switch(\n    method,\n    block    = estimate_graphon_block(A, candidate_K),\n    smooth   = estimate_graphon_smooth(A, candidate_K),\n    spectral = estimate_graphon_spectral(A, candidate_K)\n  )\n\n  out <- list(\n    W = result$W,\n    fitted = fitted_kernel_matrix(result$W, result$labels),\n    K = candidate_K,\n    n = nrow(A),\n    node_labels = result$labels,\n    method = method,\n    block_sizes = table(result$labels),\n    original_matrix = A,\n    raw_total_weight = sum(A),\n    normalization = \"none\",\n    K_selection = if (is.null(K)) \"adaptive_snapshot\" else \"fixed_or_supplied\"\n  )\n  class(out) <- \"directed_graphon\"\n  out\n}\n\nget_block_membership <- function(graphon, node_names = NULL,\n                                 output_format = c(\"list\", \"data.frame\")) {\n  if (!inherits(graphon, \"directed_graphon\")) {\n    stop(\"Input must be a directed_graphon object.\")\n  }\n  output_format <- match.arg(output_format)\n  ids <- if (is.null(node_names)) seq_len(graphon$n) else node_names\n  if (length(ids) != graphon$n) stop(\"node_names has incorrect length.\")\n\n  if (output_format == \"data.frame\") {\n    out <- data.frame(node = ids, block = graphon$node_labels, stringsAsFactors = FALSE)\n    out <- out[order(out$block, out$node), , drop = FALSE]\n    rownames(out) <- NULL\n    return(out)\n  }\n\n  out <- split(ids, factor(graphon$node_labels, levels = seq_len(graphon$K)))\n  names(out) <- paste0(\"Block_\", seq_len(graphon$K))\n  out\n}\n\ntop_singular_values <- function(A, k) {\n  A <- as_weight_matrix(A)\n  if (k < 1L) return(numeric())\n  tryCatch(\n    {\n      if (nrow(A) > 50L) {\n        as.numeric(RSpectra::svds(A, k = k)$d)\n      } else {\n        as.numeric(svd(A, nu = 0, nv = 0)$d[seq_len(k)])\n      }\n    },\n    error = function(e) as.numeric(svd(A, nu = 0, nv = 0)$d[seq_len(k)])\n  )\n}\n\nsingular_spectrum_distance <- function(A1, A2, k = 10L) {\n  A1 <- as_weight_matrix(A1)\n  A2 <- as_weight_matrix(A2)\n  if (!identical(dim(A1), dim(A2))) {\n    stop(\"Adjacency matrices must have equal dimensions.\")\n  }\n\n  k <- min(as.integer(k), nrow(A1) - 1L, nrow(A2) - 1L)\n  if (k < 1L) return(0)\n\n  s1 <- top_singular_values(A1, k)\n  s2 <- top_singular_values(A2, k)\n  sqrt(sum((s1 - s2)^2)) / sqrt(k)\n}\n\nspectral_distance <- function(A1, A2, k = 10L) {\n  singular_spectrum_distance(A1, A2, k = k)\n}\n\nas_fitted_matrix <- function(x) {\n  if (inherits(x, \"directed_graphon\")) return(x$fitted)\n  as.matrix(x)\n}\n\nfrobenius_distance <- function(M1, M2, A1 = NULL, A2 = NULL) {\n  if (is.null(M1) || is.null(M2)) {\n    if (is.null(A1) || is.null(A2)) return(NA_real_)\n    M1 <- as_weight_matrix(A1)\n    M2 <- as_weight_matrix(A2)\n  } else {\n    M1 <- as_fitted_matrix(M1)\n    M2 <- as_fitted_matrix(M2)\n  }\n  if (!identical(dim(M1), dim(M2))) stop(\"Fitted kernels must have equal dimensions.\")\n  norm(M1 - M2, type = \"F\") / nrow(M1)\n}\n\ncut_distance <- function(M1, M2, A1 = NULL, A2 = NULL, num_samples = 1000L,
                         cut_restarts = 200L) {\n  if (is.null(M1) || is.null(M2)) {\n    if (is.null(A1) || is.null(A2)) return(NA_real_)\n    M1 <- as_weight_matrix(A1)\n    M2 <- as_weight_matrix(A2)\n  } else {\n    M1 <- as_fitted_matrix(M1)\n    M2 <- as_fitted_matrix(M2)\n  }\n  if (!identical(dim(M1), dim(M2))) stop(\"Fitted kernels must have equal dimensions.\")\n\n  D <- M1 - M2\n  n <- nrow(D)\n  set.seed(42)\n\n  # FIX (2026-08): the cut norm is a maximum over all subset pairs (S, T).\n  # Drawing num_samples random pairs from a 2^n x 2^n space is a very loose\n  # lower bound: on these data it disagreed with a proper search at rank\n  # correlation 0.34 and changed which year pair ranked first.\n  #\n  # Alternating maximisation is used instead. For fixed S the optimal T is\n  # available in closed form (columns whose column sum over S is positive)\n  # and vice versa, so the iteration increases the objective monotonically\n  # and terminates. Deterministic and random restarts are run for both +D\n  # and -D. The original random-sampling result is retained as a floor, so\n  # the returned value is never below the previous implementation.\n  #\n  # NOTE: this is a tighter lower bound, not an exact value. Exact\n  # computation of the cut norm is NP-hard. Report it as\n  # \"alternating maximisation with N random restarts\".\n  ascend <- function(M, S) {\n    for (it in seq_len(100L)) {\n      cs <- colSums(M[S, , drop = FALSE])\n      Tn <- cs > 0\n      if (!any(Tn)) Tn <- cs >= max(cs)\n      rs <- rowSums(M[, Tn, drop = FALSE])\n      Sn <- rs > 0\n      if (!any(Sn)) Sn <- rs >= max(rs)\n      if (identical(Sn, S)) return(sum(M[Sn, Tn, drop = FALSE]))\n      S <- Sn\n    }\n    cs <- colSums(M[S, , drop = FALSE])\n    Tn <- cs > 0\n    if (!any(Tn)) Tn <- cs >= max(cs)\n    sum(M[S, Tn, drop = FALSE])\n  }\n\n  max_cut <- 0\n  starts_det <- list(rep(TRUE, n), rowSums(D) > 0, rowSums(D) < 0)\n  for (M in list(D, -D)) {\n    for (S0 in starts_det) {\n      if (any(S0)) max_cut <- max(max_cut, ascend(M, S0))\n    }\n    for (r in seq_len(as.integer(cut_restarts))) {\n      S0 <- sample(c(TRUE, FALSE), n, replace = TRUE)\n      if (any(S0)) max_cut <- max(max_cut, ascend(M, S0))\n    }\n  }\n\n  # Floor from the previous random-sampling estimator.\n  for (i in seq_len(as.integer(num_samples))) {\n    S <- sample(c(TRUE, FALSE), n, replace = TRUE)\n    T <- sample(c(TRUE, FALSE), n, replace = TRUE)\n    max_cut <- max(max_cut, abs(sum(D[S, T, drop = FALSE])))\n  }\n  max_cut / (n * n)\n}\n\nwasserstein_distance <- function(M1, M2) {\n  M1 <- as_fitted_matrix(M1)\n  M2 <- as_fitted_matrix(M2)\n  if (length(M1) != length(M2)) stop(\"Fitted kernels must have equal sizes.\")\n  mean(abs(sort(as.numeric(M1)) - sort(as.numeric(M2))))\n}\n\njs_divergence <- function(M1, M2) {\n  M1 <- as_fitted_matrix(M1)\n  M2 <- as_fitted_matrix(M2)\n  p <- pmax(as.numeric(M1), 0)\n  q <- pmax(as.numeric(M2), 0)\n  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)\n\n  p <- p / sum(p)\n  q <- q / sum(q)\n  m <- 0.5 * (p + q)\n  kl <- function(a, b) {\n    keep <- a > 0 & b > 0\n    sum(a[keep] * log(a[keep] / b[keep]))\n  }\n  0.5 * kl(p, m) + 0.5 * kl(q, m)\n}\n\ntv_distance <- function(M1, M2) {\n  M1 <- as_fitted_matrix(M1)\n  M2 <- as_fitted_matrix(M2)\n  p <- pmax(as.numeric(M1), 0)\n  q <- pmax(as.numeric(M2), 0)\n  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)\n\n  p <- p / sum(p)\n  q <- q / sum(q)\n  0.5 * sum(abs(p - q))\n}\n\ncompute_graphon_distance <- function(graphon1, graphon2,\n                                     method = c(\"all\", \"frobenius\", \"cut\", \"wasserstein\",\n                                                \"js\", \"tv\", \"singular_spectrum\", \"spectral\")) {\n  method <- match.arg(method)\n  F1 <- graphon1$fitted\n  F2 <- graphon2$fitted\n  A1 <- graphon1$original_matrix\n  A2 <- graphon2$original_matrix\n\n  if (!identical(dim(F1), dim(F2))) {\n    stop(\"Graphons must refer to the same node set for node-aligned comparison.\")\n  }\n\n  if (method == \"all\") {\n    return(list(\n      singular_spectrum = singular_spectrum_distance(A1, A2),\n      frobenius = frobenius_distance(F1, F2),\n      cut = cut_distance(F1, F2),\n      wasserstein = wasserstein_distance(F1, F2),\n      jensen_shannon = js_divergence(F1, F2),\n      total_variation = tv_distance(F1, F2)\n    ))\n  }\n  if (method %in% c(\"singular_spectrum\", \"spectral\")) return(singular_spectrum_distance(A1, A2))\n  if (method == \"frobenius\") return(frobenius_distance(F1, F2))\n  if (method == \"cut\") return(cut_distance(F1, F2))\n  if (method == \"wasserstein\") return(wasserstein_distance(F1, F2))\n  if (method == \"js\") return(js_divergence(F1, F2))\n  if (method == \"tv\") return(tv_distance(F1, F2))\n}\n\ngraphon_distance <- compute_graphon_distance\n\nresolve_K_schedule <- function(network_list, K = NULL) {\n  n_net <- length(network_list)\n  if (is.null(K)) {\n    out <- vapply(network_list, function(x) estimate_number_of_blocks(as_weight_matrix(x)), integer(1))\n    return(out)\n  }\n\n  K <- as.integer(K)\n  if (length(K) == 1L) return(rep.int(K, n_net))\n  if (length(K) != n_net) {\n    stop(\"K must be NULL, one integer, or one integer per network.\")\n  }\n  K\n}\n\ncompare_graphons <- function(network_list, K = NULL,\n                             method = c(\"block\", \"smooth\", \"spectral\")) {\n  method <- match.arg(method)\n  if (length(network_list) < 2L) stop(\"At least two networks are required.\")\n\n  network_names <- names(network_list)\n  if (is.null(network_names)) network_names <- as.character(seq_along(network_list))\n\n  dims <- vapply(network_list, function(x) nrow(as_weight_matrix(x)), integer(1))\n  if (length(unique(dims)) != 1L) {\n    stop(\"All networks in a comparison set must contain the same number of nodes.\")\n  }\n\n  K_schedule <- resolve_K_schedule(network_list, K = K)\n  names(K_schedule) <- network_names\n\n  graphons <- lapply(seq_along(network_list), function(i) {\n    estimate_directed_graphon(network_list[[i]], K = K_schedule[i], method = method)\n  })\n  names(graphons) <- network_names\n\n  metric_names <- c(\n    \"singular_spectrum\", \"frobenius\", \"cut\", \"wasserstein\",\n    \"jensen_shannon\", \"total_variation\"\n  )\n  distances <- setNames(lapply(metric_names, function(x) {\n    matrix(0, length(graphons), length(graphons),\n           dimnames = list(network_names, network_names))\n  }), metric_names)\n\n  for (i in seq_len(length(graphons) - 1L)) {\n    for (j in seq.int(i + 1L, length(graphons))) {\n      d <- compute_graphon_distance(graphons[[i]], graphons[[j]], method = \"all\")\n      for (nm in metric_names) {\n        distances[[nm]][i, j] <- distances[[nm]][j, i] <- d[[nm]]\n      }\n    }\n  }\n\n  distances$spectral <- distances$singular_spectrum\n\n  list(\n    graphons = graphons,\n    distances = distances,\n    network_names = network_names,\n    K = K_schedule,\n    K_by_network = K_schedule,\n    method = method,\n    normalization = \"none\",\n    distance_domain = \"node_aligned_fitted_kernel\",\n    benchmark = \"singular_spectrum_on_original_weighted_adjacency\"\n  )\n}\n\nplot_graphon <- function(graphon, title = NULL) {\n  if (!inherits(graphon, \"directed_graphon\")) {\n    stop(\"graphon must be a directed_graphon object.\")\n  }\n  df <- reshape2::melt(graphon$W)\n  names(df) <- c(\"from_block\", \"to_block\", \"weight\")\n  ggplot2::ggplot(df, ggplot2::aes(x = to_block, y = from_block, fill = weight)) +\n    ggplot2::geom_tile() +\n    ggplot2::scale_y_reverse() +\n    ggplot2::labs(\n      x = \"Destination block\", y = \"Origin block\", fill = \"Mean weight\",\n      title = if (is.null(title)) paste0(\"Directed weighted block kernel (K=\", graphon$K, \")\") else title\n    ) +\n    ggplot2::theme_minimal()\n}\n\nplot_distance_heatmap <- function(result, metric = \"frobenius\") {\n  D <- result$distances[[metric]]\n  if (is.null(D)) stop(\"Unknown distance metric: \", metric)\n  df <- reshape2::melt(D)\n  names(df) <- c(\"network_1\", \"network_2\", \"distance\")\n  ggplot2::ggplot(df, ggplot2::aes(x = network_2, y = network_1, fill = distance)) +\n    ggplot2::geom_tile() +\n    ggplot2::labs(x = NULL, y = NULL, fill = metric) +\n    ggplot2::theme_minimal() +\n    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))\n}\n\ngenerate_pastel_palette <- function(n) {\n  if (n <= 0L) return(character())\n  hues <- seq(15, 375, length.out = n + 1L)[seq_len(n)]\n  grDevices::hcl(h = hues, c = 45, l = 75)\n}\n\nmix_colors_additive <- function(colors) {\n  if (length(colors) == 0L) return(\"gray70\")\n  rgb_mtx <- grDevices::col2rgb(colors) / 255\n  m <- pmin(1, rowMeans(rgb_mtx))\n  grDevices::rgb(m[1], m[2], m[3])\n}\n\ndarken_color <- function(color, factor = 0.75) {\n  rgb <- grDevices::col2rgb(color) * factor\n  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)\n}\n\nannotate_graph_with_graphon <- function(g,\n                                        K = NULL,\n                                        method = \"block\",\n                                        block_colors = NULL,\n                                        default_color = \"gray70\",\n                                        log_transform = c(\"none\", \"log\", \"loglog\"),\n                                        min_width = 0.5,\n                                        max_width = 6) {\n  log_transform <- match.arg(log_transform)\n  if (!inherits(g, \"igraph\")) stop(\"g must be an igraph object.\")\n\n  A <- as_weight_matrix(g)\n  graphon_obj <- estimate_directed_graphon(A, K = K, method = method)\n  member <- graphon_obj$node_labels\n\n  g <- igraph::set_vertex_attr(g, \"graphonmember\", value = member)\n  out_strength <- rowSums(A)\n  in_strength <- colSums(A)\n  g <- igraph::set_vertex_attr(\n    g, \"origmember\", value = ifelse(out_strength > 0, member, NA_integer_)\n  )\n  g <- igraph::set_vertex_attr(\n    g, \"destmember\", value = ifelse(in_strength > 0, member, NA_integer_)\n  )\n\n  if (is.null(block_colors)) block_colors <- generate_pastel_palette(graphon_obj$K)\n  if (length(block_colors) < graphon_obj$K) {\n    block_colors <- rep(block_colors, length.out = graphon_obj$K)\n  }\n  if (length(block_colors) > 0L) {\n    g <- igraph::set_vertex_attr(g, \"color\", value = block_colors[member])\n  } else {\n    g <- igraph::set_vertex_attr(g, \"color\", value = rep(default_color, igraph::vcount(g)))\n  }\n\n  if (igraph::ecount(g) > 0L) {\n    ends <- igraph::ends(g, igraph::E(g), names = FALSE)\n    edge_cols <- vapply(seq_len(nrow(ends)), function(i) {\n      mix_colors_additive(block_colors[member[ends[i, ]]])\n    }, character(1))\n    g <- igraph::set_edge_attr(g, \"color\", value = edge_cols)\n\n    w <- if (\"weight\" %in% igraph::edge_attr_names(g)) {\n      as.numeric(igraph::edge_attr(g, \"weight\"))\n    } else {\n      rep(1, igraph::ecount(g))\n    }\n    if (log_transform == \"log\") w <- log1p(pmax(w, 0))\n    if (log_transform == \"loglog\") w <- log1p(log1p(pmax(w, 0)))\n\n    finite_w <- w[is.finite(w)]\n    if (length(finite_w) == 0L || length(unique(finite_w)) <= 1L) {\n      widths <- rep((min_width + max_width) / 2, length(w))\n    } else {\n      r <- range(finite_w)\n      widths <- min_width + (w - r[1]) / (r[2] - r[1]) * (max_width - min_width)\n      widths[!is.finite(widths)] <- min_width\n    }\n    g <- igraph::set_edge_attr(g, \"width\", value = widths)\n  }\n\n  attr(g, \"graphon\") <- graphon_obj\n  g\n}\n\n"
GEOPLOT_HELPER_SOURCE <- "geoplot <- function(g,\n                    mode = c(\"basemap\", \"rworldmap\", \"custom_shape\"),\n                    basemap = \"OpenStreetMap\",\n                    shape = NULL,\n                    show_tmap_on_custom = FALSE,\n                    node_col = \"black\",\n                    node_alpha = 0.8,\n                    node_size = 0.8,\n                    node_size_scale = c(0.5, 2.5),\n                    colormap = grDevices::rainbow,\n                    edge_col = \"gray\",\n                    edge_alpha = 0.5,\n                    edge_width = 1,\n                    edge_curved = FALSE,\n                    curvature = 0.15,\n                    arrow_size = 0.0,\n                    shape_bg = \"#f0f0f0\",\n                    shape_border_col = \"#d0d0d0\",\n                    # [NEW] Label parameters\n                    show_labels = FALSE,\n                    label_field = \"nodeLabel\",\n                    label_col = \"black\",\n                    label_size_mult = 1,\n                    top_n_labels = NULL) {\n\n  library(igraph); library(sf); library(tmap); library(scales)\n  #tmap_options(facet.max = 2000)\n  mode <- match.arg(mode)\n\n  v_count <- igraph::vcount(g)\n  e_count <- igraph::ecount(g)\n  nodes_df <- igraph::as_data_frame(g, what = \"vertices\")\n\n  # [NEW] Detect whether the graph is directed (controls arrowhead rendering)\n  is_dir <- igraph::is_directed(g)\n\n  # Detect pre-existing edge attributes on the igraph object\n  has_edge_color_attr <- e_count > 0 && !is.null(igraph::edge_attr(g, \"color\"))\n  has_edge_width_attr <- e_count > 0 && !is.null(igraph::edge_attr(g, \"width\"))\n\n  # --- 1. NODE DATA (color: handles functions and fixed palettes) ---\n  get_colors <- function(n_levels) {\n    if (is.function(colormap)) {\n      return(colormap(n_levels))\n    } else if (is.character(colormap) && length(colormap) >= n_levels) {\n      return(colormap[1:n_levels])\n    } else {\n      return(grDevices::rainbow(n_levels))\n    }\n  }\n\n  if (length(node_col) == 1 && node_col %in% names(nodes_df)) {\n    vals <- as.factor(nodes_df[[node_col]])\n    node_colors_raw <- get_colors(length(levels(vals)))[vals]\n  } else if (length(node_col) == v_count && !is.character(node_col)) {\n    vals <- as.factor(node_col)\n    node_colors_raw <- get_colors(length(levels(vals)))[vals]\n  } else if (length(node_col) == v_count) {\n    node_colors_raw <- node_col\n  } else {\n    node_colors_raw <- rep(node_col[1], v_count)\n  }\n  final_node_colors <- scales::alpha(node_colors_raw, node_alpha)\n\n  # --- NODE SIZE ---\n  rescale_vec <- function(x, range_to) {\n    if (length(x) == 0) return(numeric(0))\n    if (all(x == x[1])) return(rep(range_to[1], length(x)))\n    x_num <- as.numeric(x)\n    (x_num - min(x_num, na.rm = TRUE)) /\n      (max(x_num, na.rm = TRUE) - min(x_num, na.rm = TRUE)) *\n      (range_to[2] - range_to[1]) + range_to[1]\n  }\n\n  # [CHANGED] Guard: coerce list to atomic vector (igraph may store attrs as lists)\n  if (is.list(node_size) && !is.data.frame(node_size)) {\n    node_size <- unlist(node_size)\n  }\n\n  if (is.character(node_size) && length(node_size) == 1 && node_size %in% names(nodes_df)) {\n    col_data <- nodes_df[[node_size]]\n    if (is.list(col_data)) col_data <- unlist(col_data)\n    final_node_size <- rescale_vec(as.numeric(col_data), node_size_scale)\n  } else if (length(node_size) == v_count) {\n    final_node_size <- rescale_vec(as.numeric(node_size), node_size_scale)\n  } else {\n    size_val <- suppressWarnings(as.numeric(node_size[1]))\n    if (is.na(size_val)) {\n      warning(\"node_size could not be interpreted as numeric. Falling back to 0.8.\")\n      size_val <- 0.8\n    }\n    final_node_size <- rep(size_val, v_count)\n  }\n\n  nodes_sf <- sf::st_as_sf(nodes_df, coords = c(\"lon\", \"lat\"), crs = 4326)\n  nodes_sf$node_color_val <- final_node_colors\n  nodes_sf$node_size_val  <- final_node_size * 0.05\n\n  # --- [NEW] NODE LABELS ---\n  # Precompute label text, font sizes, and the display mask so that both\n  # the tmap and the static rendering paths can use them.\n  label_text <- NULL\n  label_cex  <- NULL\n  label_mask <- NULL\n\n  if (show_labels) {\n    # Resolve the label text source\n    if (label_field %in% names(nodes_df)) {\n      label_text <- as.character(nodes_df[[label_field]])\n    } else if (\"name\" %in% names(nodes_df)) {\n      warning(paste0(\"Label field '\", label_field, \"' not found in vertex \",\n                     \"attributes. Falling back to 'name'.\"))\n      label_text <- as.character(nodes_df[[\"name\"]])\n    } else {\n      warning(paste0(\"Label field '\", label_field, \"' not found and no 'name' \",\n                     \"attribute available. Using node indices.\"))\n      label_text <- as.character(seq_len(v_count))\n    }\n\n    # Font size proportional to node size, scaled by user multiplier\n    label_cex <- final_node_size * label_size_mult\n\n    # Optionally restrict to the top-N largest nodes\n    if (!is.null(top_n_labels) && is.numeric(top_n_labels) &&\n        top_n_labels > 0 && top_n_labels < v_count) {\n      top_idx    <- order(final_node_size, decreasing = TRUE)[1:top_n_labels]\n      label_mask <- seq_len(v_count) %in% top_idx\n    } else {\n      label_mask <- rep(TRUE, v_count)\n    }\n\n    # Store in nodes_sf for tmap modes\n    nodes_sf$node_label    <- ifelse(label_mask, label_text, NA_character_)\n    nodes_sf$label_size_val <- ifelse(label_mask, label_cex * 0.4, 0)\n  }\n\n  # --- 2. EDGE DATA ---\n  edges_sf  <- NULL\n  arrows_sf <- NULL          # [NEW] Arrowhead polygons for directed graphs\n\n  if (e_count > 0) {\n    edge_list <- igraph::as_edgelist(g, names = TRUE)\n\n    # -- [CHANGED] Edge colors: prefer E(g)$color when available --\n    if (has_edge_color_attr) {\n      e_cols_raw <- igraph::E(g)$color\n      if (is.list(e_cols_raw)) e_cols_raw <- unlist(e_cols_raw)\n    } else if (edge_col == \"mixed\") {\n      from_idx   <- match(edge_list[, 1], nodes_df$name)\n      e_cols_raw <- node_colors_raw[from_idx]\n    } else {\n      e_cols_raw <- rep(edge_col[1], e_count)\n    }\n    final_edge_colors <- scales::alpha(e_cols_raw, edge_alpha)\n\n    # -- [NEW] Edge widths: prefer E(g)$width when available --\n    if (has_edge_width_attr) {\n      final_edge_widths <- igraph::E(g)$width\n      if (is.list(final_edge_widths)) final_edge_widths <- unlist(final_edge_widths)\n      final_edge_widths <- as.numeric(final_edge_widths)\n      final_edge_widths[is.na(final_edge_widths)] <- edge_width\n    } else {\n      final_edge_widths <- rep(edge_width, e_count)\n    }\n\n    # -- Edge geometries (straight or Bezier) --\n    calculate_bezier <- function(p1, p2, v, n = 30) {\n      mid <- (p1 + p2) / 2\n      diff <- p2 - p1\n      norm <- c(-diff[2], diff[1])\n      control <- mid + norm * v\n      t <- seq(0, 1, length.out = n)\n      coords <- matrix(NA, n, 2)\n      for (i in 1:n)\n        coords[i, ] <- (1 - t[i])^2 * p1 +\n        2 * (1 - t[i]) * t[i] * control +\n        t[i]^2 * p2\n      return(coords)\n    }\n\n    edge_geoms <- lapply(1:e_count, function(i) {\n      p1 <- as.numeric(nodes_df[nodes_df$name == edge_list[i, 1], c(\"lon\", \"lat\")])\n      p2 <- as.numeric(nodes_df[nodes_df$name == edge_list[i, 2], c(\"lon\", \"lat\")])\n      sf::st_linestring(\n        if (edge_curved) calculate_bezier(p1, p2, curvature) else rbind(p1, p2)\n      )\n    })\n\n    # [CHANGED] edges_sf now carries per-edge width as well\n    edges_sf <- sf::st_sf(\n      geometry       = sf::st_sfc(edge_geoms, crs = 4326),\n      edge_color_val = as.character(final_edge_colors),\n      edge_width_val = final_edge_widths\n    )\n\n    # ------------------------------------------------------------------\n    # [NEW] Build arrowhead polygons for directed graphs when arrow_size > 0\n    # ------------------------------------------------------------------\n    # Arrowheads are drawn as filled triangle polygons at the target end of\n    # each edge.  The size is relative to the spatial extent of the node\n    # bounding box so that arrowheads look proportional regardless of zoom.\n    # arrow_size acts as a user-controlled multiplier (default 0 = no arrows).\n    # ------------------------------------------------------------------\n    if (is_dir && arrow_size > 0) {\n\n      # Compute a base arrowhead length from the spatial extent of all nodes\n      bbox <- sf::st_bbox(nodes_sf)\n      extent_diag <- sqrt((bbox[\"xmax\"] - bbox[\"xmin\"])^2 +\n                            (bbox[\"ymax\"] - bbox[\"ymin\"])^2)\n      # Fallback for degenerate cases (all nodes at the same location)\n      if (!is.finite(extent_diag) || extent_diag == 0) extent_diag <- 1\n\n      # Base length = 1.5 % of the diagonal, scaled by the user's arrow_size\n      base_arrow_len <- extent_diag * 0.015 * arrow_size\n\n      # Helper: create a single arrowhead triangle polygon\n      # tip       - numeric(2), the point of the arrowhead (target end)\n      # direction - numeric(2), vector pointing from penultimate to tip\n      # size      - scalar, length of the arrowhead along the direction axis\n      make_arrowhead <- function(tip, direction, size) {\n        d_len <- sqrt(sum(direction^2))\n        if (d_len == 0) return(NULL)\n        d    <- direction / d_len          # unit vector along edge\n        perp <- c(-d[2], d[1])             # perpendicular unit vector\n        base_mid <- tip - d * size         # midpoint of the triangle base\n        p1 <- base_mid + perp * size * 0.45\n        p2 <- base_mid - perp * size * 0.45\n        sf::st_polygon(list(rbind(tip, p1, p2, tip)))\n      }\n\n      arrow_polys  <- vector(\"list\", e_count)\n      arrow_colors <- character(e_count)\n      keep         <- logical(e_count)\n\n      for (i in seq_len(e_count)) {\n        coords <- sf::st_coordinates(edges_sf$geometry[i])\n        n_pts  <- nrow(coords)\n        # Tip = last point (target node), prev = penultimate point\n        tip  <- coords[n_pts, 1:2]\n        prev <- coords[max(n_pts - 1, 1), 1:2]\n        direction <- tip - prev\n\n        poly <- make_arrowhead(as.numeric(tip),\n                               as.numeric(direction),\n                               base_arrow_len)\n        if (!is.null(poly)) {\n          arrow_polys[[i]] <- poly\n          arrow_colors[i]  <- edges_sf$edge_color_val[i]\n          keep[i]          <- TRUE\n        }\n      }\n\n      if (any(keep)) {\n        arrows_sf <- sf::st_sf(\n          geometry    = sf::st_sfc(arrow_polys[keep], crs = 4326),\n          arrow_color = arrow_colors[keep]\n        )\n      }\n    }\n  }\n\n  # --- 3. RENDERING ---\n  if (mode == \"basemap\" || (mode == \"custom_shape\" && show_tmap_on_custom)) {\n    # ---- Interactive (tmap) mode ----\n    tmap_mode(\"view\")\n\n    final_bmap <- if (basemap == \"osm\") \"OpenStreetMap\"\n    else if (basemap == \"dark\") \"CartoDB.DarkMatter\"\n    else if (basemap == \"light\") \"CartoDB.Positron\"\n    else basemap\n\n    mapa <- if (mode == \"basemap\") {\n      tm_basemap(final_bmap)\n    } else {\n      tm_shape(st_transform(shape, 4326)) +\n        tm_polygons(col = shape_bg, border.col = shape_border_col)\n    }\n\n    # [CHANGED] Edges: variable width from edge_width_val column\n    if (!is.null(edges_sf)) {\n      mapa <- mapa +\n        tm_shape(edges_sf) +\n        tm_lines(col     = \"edge_color_val\",\n                 lwd     = \"edge_width_val\",\n                 col.legend = tm_legend_hide(),\n                 lwd.legend = tm_legend_hide())\n    }\n\n    # [NEW] Arrowheads (tmap): rendered as filled triangle polygons on top\n    # of the edge lines so that the direction of each edge is visible.\n    if (!is.null(arrows_sf)) {\n      mapa <- mapa +\n        tm_shape(arrows_sf) +\n        tm_polygons(col        = \"arrow_color\",\n                    border.col = \"arrow_color\",\n                    col.legend = tm_legend_hide())\n    }\n\n    # Nodes\n    mapa <- mapa +\n      tm_shape(nodes_sf) +\n      tm_symbols(col    = \"node_color_val\",\n                 size   = \"node_size_val\",\n                 border.lwd = 0,\n                 col.legend  = tm_legend_hide(),\n                 size.legend = tm_legend_hide())\n\n    # [NEW] Labels (tmap)\n    if (show_labels) {\n      # Subset to labeled nodes only (non-NA label)\n      labeled_sf <- nodes_sf[!is.na(nodes_sf$node_label), ]\n      if (nrow(labeled_sf) > 0) {\n        mapa <- mapa +\n          tm_shape(labeled_sf) +\n          tm_text(\"node_label\",\n                  size = \"label_size_val\",\n                  col  = label_col,\n                  size.legend = tm_legend_hide())\n      }\n    }\n\n    return(mapa)\n\n  } else {\n    # ---- Static (base R) mode ----\n    if (mode == \"rworldmap\") {\n      world     <- rworldmap::getMap(resolution = \"high\")\n      shape_gps <- sf::st_as_sf(world)\n      bbox      <- st_bbox(nodes_sf)\n      lims      <- list(x = c(bbox[\"xmin\"], bbox[\"xmax\"]),\n                        y = c(bbox[\"ymin\"], bbox[\"ymax\"]))\n    } else {\n      shape_gps <- st_transform(shape, 4326)\n      lims      <- NULL\n    }\n\n    plot(st_geometry(shape_gps),\n         col = shape_bg, border = shape_border_col,\n         xlim = lims$x, ylim = lims$y)\n\n    # [CHANGED] Edges: per-edge color AND per-edge width\n    if (!is.null(edges_sf)) {\n      plot(st_geometry(edges_sf), add = TRUE,\n           col = edges_sf$edge_color_val,\n           lwd = edges_sf$edge_width_val)\n    }\n\n    # [NEW] Arrowheads (base R): filled triangle polygons overlaid on edges\n    if (!is.null(arrows_sf)) {\n      plot(st_geometry(arrows_sf), add = TRUE,\n           col    = arrows_sf$arrow_color,\n           border = arrows_sf$arrow_color)\n    }\n\n    # Nodes\n    plot(st_geometry(nodes_sf), add = TRUE,\n         col = nodes_sf$node_color_val,\n         pch = 16, cex = final_node_size)\n\n    # [NEW] Labels (base R)\n    if (show_labels) {\n      idx_to_label <- which(label_mask)\n      if (length(idx_to_label) > 0) {\n        coords <- sf::st_coordinates(nodes_sf)\n        text(coords[idx_to_label, 1],\n             coords[idx_to_label, 2],\n             labels = label_text[idx_to_label],\n             cex    = label_cex[idx_to_label],\n             col    = label_col,\n             pos    = 3,\n             offset = 0.5)\n      }\n    }\n\n    return(recordPlot())\n  }\n}\n\n"
TOKABLE_HELPER_SOURCE <- "tokable <- function(D,subtitle,escape=FALSE,fontsize=8) {\n  tokable<-kable(D,digits = Inf,format=\"latex\",booktabs=TRUE,longtable=FALSE,escape = escape,linesep=\"\") %>% row_spec(0,bold=TRUE) %>% column_spec(1,bold=TRUE)  %>% footnote(general=paste(\"\\\\\\\\\\\\centering\",subtitle,sep=\" \"),general_title = \"\", threeparttable = TRUE,escape = escape) %>% kable_styling(latex_options = c(\"striped\",\"repeat_header\"),full_width = FALSE,font_size = fontsize) %>% kable_paper()  \n}\n"


write_embedded_helper("graphon_distance_directed.R", GRAPHON_HELPER_SOURCE)
write_embedded_helper("geoplot.R", GEOPLOT_HELPER_SOURCE)
write_embedded_helper("tokable.R", TOKABLE_HELPER_SOURCE)

if (!file.exists("hungary.RData")) {
  hungary <- if (requireNamespace("maps", quietly = TRUE)) {
    tryCatch(maps::map("world", regions = "Hungary", plot = FALSE, fill = TRUE),
             error = function(e) NULL)
  } else NULL
  save(hungary, file = "hungary.RData")
  log_step("Created hungary.RData compatibility object.")
}

# =============================================================================
# Inlined metric and graphon implementations
# =============================================================================

# ---- Inlined from /TS2VG/tsnda/R/percolate.R ----
#-----------------------------------------------------------------------------#
#                                                                             #
#              Visibility Graph Based Forecast                                #
#                                                                             #
#  Written by: Zsolt T. Kosztyan                                              #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: August 2023                                                  #
#-----------------------------------------------------------------------------#

#' @export

percolate = function(g, size, d) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  giant = vector()

  # initial size of giant component
  c = igraph::components(g)
  giant[1] = max(c$csize)

  names(d) = 1:length(d)
  d = sort(d, decreasing=TRUE)
  vital = as.integer(names(d[1:size]))

  for (i in 1:size) {
    c = igraph::components(igraph::delete_vertices(g, vital[1:i]))
    giant[i+1] = max(c$csize)
  }

  return(giant)

}



# ---- Inlined from /Mob/lacunarity.R ----
packages <- c("igraph", "ggplot2", "dplyr", "gridExtra", "tidyr", "e1071")

# Install packages if needed
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set the seed for reproducibility
set.seed(123)

local_lacunarity <- function(graph, radius = 2) {
  
  n_vertices <- vcount(graph)
  
  if(n_vertices <= 2) {
    return(rep(1, n_vertices))
  }
  
  dist_matrix <- distances(graph)
  local_lac <- numeric(n_vertices)
  
  for(i in 1:n_vertices) {
    # Neighborhood of node i
    neighbors_in_radius <- which(dist_matrix[i, ] <= radius &
                                   dist_matrix[i, ] < Inf)
    
    if(length(neighbors_in_radius) > 1) {
      # Degrees of neighboring nodes (local masses)
      neighbor_degrees <- igraph::degree(graph)[neighbors_in_radius]
      
      if(var(neighbor_degrees) > 0) {
        mean_deg <- mean(neighbor_degrees)
        mean_deg_squared <- mean(neighbor_degrees^2)
        local_lac[i] <- mean_deg_squared / (mean_deg^2)
      } else {
        local_lac[i] <- 1
      }
    } else {
      local_lac[i] <- 1
    }
  }
  
  return(local_lac)
}


network_lacunarity <- function(graph, radii = NULL, n_samples = 20) {
  
  n_vertices <- vcount(graph)
  
  # Return a trivial result when the graph is too small
  if(n_vertices <= 2) {
    return(list(lacunarity = rep(1, length(radii %||% c(1))),
                radii = radii %||% c(1), masses = list(c(1))))
  }
  
  # Compute the distance matrix
  dist_matrix <- distances(graph)
  max_distance <- max(dist_matrix[is.finite(dist_matrix)])
  
  if(is.null(radii)) {
    if(max_distance == 0) max_distance <- 2
    radii <- 1:min(max_distance, 8)
  }
  
  # Select sampling centers
  if(n_samples > n_vertices) n_samples <- n_vertices
  sample_centers <- sample(1:n_vertices, min(n_samples, n_vertices))
  
  lacunarity_values <- numeric(length(radii))
  all_masses <- list()
  
  for(i in seq_along(radii)) {
    r <- radii[i]
    masses <- numeric(length(sample_centers))
    
    # Compute the mass for every sampling center
    for(j in seq_along(sample_centers)) {
      center <- sample_centers[j]
      distances_from_center <- dist_matrix[center, ]
      
      # Number of nodes within distance r (mass)
      masses[j] <- sum(distances_from_center <= r & is.finite(distances_from_center))
    }
    
    all_masses[[i]] <- masses
    
    # Compute lacunarity: <M^2> / <M>^2
    if(length(masses) > 1 && var(masses) > 0) {
      mean_mass <- mean(masses)
      mean_mass_squared <- mean(masses^2)
      
      if(mean_mass > 0) {
        lacunarity_values[i] <- mean_mass_squared / (mean_mass^2)
      } else {
        lacunarity_values[i] <- 1
      }
    } else {
      lacunarity_values[i] <- 1
    }
  }
  
  return(list(
    lacunarity = lacunarity_values,
    radii = radii,
    masses = all_masses,
    sample_size = length(sample_centers)
  ))
}

network_box_counting <- function(graph, box_sizes = NULL) {
  
  n_vertices <- vcount(graph)
  
  # Return a trivial result when the graph is empty or too small
  if(n_vertices <= 1) {
    return(list(dimension = NA, box_sizes = c(1), n_boxes = c(1)))
  }
  
  # Compute the distance matrix
  dist_matrix <- distances(graph)
  diameter_graph <- max(dist_matrix[is.finite(dist_matrix)])
  
  if(is.null(box_sizes)) {
    if(diameter_graph == 0) diameter_graph <- 1
    box_sizes <- seq(1, max(3, diameter_graph), by = 1)
  }
  
  n_boxes <- numeric(length(box_sizes))
  
  for(i in seq_along(box_sizes)) {
    l_B <- box_sizes[i]  # Box size
    
    # Group nodes according to box_size
    uncovered <- 1:n_vertices
    boxes <- 0
    
    while(length(uncovered) > 0) {
      # Select a center node
      center <- uncovered[1]
      
      # Find all nodes within distance l_B
      distances_from_center <- dist_matrix[center, uncovered]
      in_box <- which(distances_from_center <= l_B)
      
      # Remove these nodes
      uncovered <- uncovered[-in_box]
      boxes <- boxes + 1
    }
    
    n_boxes[i] <- boxes
  }
  
  # Log-log regression for estimating the dimension
  # N(l_B) ~ l_B^(-d_B), where d_B is the box-counting dimension
  valid_idx <- n_boxes > 0 & box_sizes > 0 & is.finite(n_boxes)
  
  if(sum(valid_idx) < 3) {
    return(list(dimension = NA, box_sizes = box_sizes, n_boxes = n_boxes, r_squared = NA))
  }
  
  log_box_sizes <- log(box_sizes[valid_idx])
  log_n_boxes <- log(n_boxes[valid_idx])
  
  # Check whether there is any variation
  if(var(log_box_sizes) == 0 || var(log_n_boxes) == 0) {
    return(list(dimension = NA, box_sizes = box_sizes, n_boxes = n_boxes, r_squared = NA))
  }
  
  fit <- lm(log_n_boxes ~ log_box_sizes)
  dimension <- -as.numeric(coef(fit)[2])  # Negative because N(l_B) ~ l_B^(-d_B)
  
  return(list(
    dimension = dimension,
    box_sizes = box_sizes,
    n_boxes = n_boxes,
    fit = fit,
    r_squared = summary(fit)$r.squared
  ))
}

sandbox_dimension <- function(graph, radii = NULL, n_centers = 10) {
  
  n_vertices <- vcount(graph)
  
  if(n_vertices <= 1) {
    return(list(dimension = NA, radii = c(1), masses = c(1), r_squared = NA))
  }
  
  # Distance matrix
  dist_matrix <- distances(graph)
  max_distance <- max(dist_matrix[is.finite(dist_matrix)])
  
  if(is.null(radii)) {
    if(max_distance == 0) max_distance <- 1
    radii <- 1:min(max_distance, 8)
  }
  
  # Select random centers
  if(n_centers > n_vertices) n_centers <- n_vertices
  centers <- sample(1:n_vertices, min(n_centers, n_vertices))
  
  masses <- matrix(0, nrow = length(centers), ncol = length(radii))
  
  for(i in seq_along(centers)) {
    center <- centers[i]
    for(j in seq_along(radii)) {
      r <- radii[j]
      # Number of nodes within distance r
      distances_from_center <- dist_matrix[center, ]
      masses[i, j] <- sum(distances_from_center <= r & is.finite(distances_from_center))
    }
  }
  
  # Mean mass for each r
  avg_masses <- colMeans(masses)
  
  # Log-log regression: M(r) ~ r^d_s
  valid_idx <- avg_masses > 0 & radii > 0 & is.finite(avg_masses)
  
  if(sum(valid_idx) < 3) {
    return(list(dimension = NA, radii = radii, masses = avg_masses, r_squared = NA))
  }
  
  log_radii <- log(radii[valid_idx])
  log_masses <- log(avg_masses[valid_idx])
  
  # Check whether there is any variation
  if(var(log_radii) == 0 || var(log_masses) == 0) {
    return(list(dimension = NA, radii = radii, masses = avg_masses, r_squared = NA))
  }
  
  fit <- lm(log_masses ~ log_radii)
  dimension <- as.numeric(coef(fit)[2])
  
  return(list(
    dimension = dimension,
    radii = radii,
    masses = avg_masses,
    fit = fit,
    r_squared = summary(fit)$r.squared
  ))
}



# ---- Inlined from /Mob/resilience_centrality.R ----
# Helper: always return a vector of labels with length = vcount(g)
safe_vertex_labels <- function(g) {
  if (!is.null(V(g)$name) && length(V(g)$name) == vcount(g)) {
    V(g)$name
  } else {
    as.character(seq_len(vcount(g)))  # fallback: 1,2,...,n
  }
}

beta_eff <- function(g) {
  deg <- igraph::degree(g, mode = "all")
  if (length(deg) == 0 || mean(deg) == 0) return(0)
  mean(deg^2) / mean(deg)
}

resilience_centrality <- function(g) {
  n <- vcount(g)
  if (n == 0) return(data.frame(vertex = character(), resilience_centrality = numeric()))
  
  base_beta <- beta_eff(g)
  crit_values <- numeric(n)
  
  for (i in seq_len(n)) {
    g_minus_i <- delete_vertices(g, i)
    beta_minus <- beta_eff(g_minus_i)
    if (base_beta == 0) {
      crit_values[i] <- 0
    } else {
      crit_values[i] <- (base_beta - beta_minus) / base_beta
    }
  }
  
  data.frame(
    vertex = safe_vertex_labels(g),
    resilience_centrality = crit_values
  )
}


# ---- Inlined from /Mob/kirchhoff_index.R ----
library(igraph)
library(Matrix)
library(RSpectra)

kirchhoff_index <- function(g, tol = 1e-12, use_symmetrized = TRUE,
                            dense_threshold = 500, max_sparse_k = 800) {
  
  # -- 1. Construct the symmetric Laplacian matrix ---------------------
  if (use_symmetrized) {
    gu <- as.undirected(g, mode = "collapse")
    L  <- laplacian_matrix(gu, sparse = TRUE)
  } else {
    Ld <- laplacian_matrix(g, sparse = TRUE)
    L  <- (Ld + t(Ld)) / 2
  }
  
  n <- nrow(L)
  if (n <= 1) return(0)
  
  # -- 2. Small and medium graphs: reliable dense eigendecomposition ---
  
  if (n <= dense_threshold) {
    evals <- eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values
    evals <- evals[evals > tol]
    return(n * sum(1 / evals))
  }
  
  # -- 3. Large graphs: shift-invert mode with RSpectra ----------------
  #
  #   which = "LM" with sigma finds the largest eigenvalues of
  #   (L - sigma I)^-1, which returns the eigenvalues closest to sigma.
  #
  #   A slightly negative sigma makes (L - sigma I) positive definite and
  #   avoids singularity. The zero eigenvalue is returned but filtered out.
  
  k     <- min(n - 2, max_sparse_k)   # n-2 leaves room for ncv
  sigma <- -1e-6                       # Small negative shift
  
  # ncv: number of ARPACK Lanczos vectors (k < ncv <= n)
  ncv <- min(max(2 * k + 1, k + 30), n)
  
  result <- tryCatch({
    
    out <- eigs_sym(L, k = k, which = "LM", sigma = sigma,
                    opts = list(ncv = ncv, maxitr = 500, tol = 1e-10))
    
    lam <- Re(out$values)
    lam <- lam[lam > tol]
    
    if (length(lam) == 0) {
      warning("No positive eigenvalues were found.")
      return(0)
    }
    
    # If k < n-1, this is a lower bound because smaller eigenvalues dominate
    # the sum of reciprocal eigenvalues.
    if (length(lam) < n - 1) {
      message(sprintf("Partial sum: %d of %d nonzero eigenvalues.",
                      length(lam), n - 1))
    }
    
    n * sum(1 / lam)
    
  }, error = function(e1) {
    
    # -- Fallback A: adjust ncv and maxitr --
    tryCatch({
      message("The first eigs_sym attempt failed; retrying with adjusted parameters...")
      ncv2 <- min(max(3 * k + 1, k + 50), n)
      k2   <- max(min(k, n - 10), 1)
      
      out <- eigs_sym(L, k = k2, which = "LM", sigma = sigma,
                      opts = list(ncv = ncv2, maxitr = 1000, tol = 1e-8))
      lam <- Re(out$values)
      lam <- lam[lam > tol]
      n * sum(1 / lam)
      
    }, error = function(e2) {
      
      # -- Fallback B: dense eigendecomposition --
      warning(paste0("RSpectra failed after retrying: ", e2$message,
                     "; using a dense eigendecomposition fallback."))
      evals <- eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values
      evals <- evals[evals > tol]
      n * sum(1 / evals)
    })
  })
  
  result
}



# ---- Inlined from /Mob/nri.R ----
library(igraph)

# Resilience Index: average fraction of nodes in LCC after random removals
resilience_index <- function(g, frac_remove = 0.2, trials = 100) {
  n <- vcount(g)
  sizes <- numeric(trials)
  
  for (i in seq_len(trials)) {
    remove <- sample(V(g), size = ceiling(frac_remove * n))
    g_minus <- delete_vertices(g, remove)
    comps <- components(g_minus)
    sizes[i] <- max(comps$csize) / (n - length(remove))
  }
  
  mean(sizes)
}

# Normalized Resilience Index: compare to complete graph baseline
normalized_resilience_index <- function(g, frac_remove = 0.2, trials = 100) {
  RI_g <- resilience_index(g, frac_remove, trials)
  # Complete graph baseline (always perfectly connected until all nodes removed)
  RI_max <- 1
  RI_g / RI_max
}



# ---- Inlined from /Mob/network_entropy.R ----
network_entropy <- function(g, mode = "all") {
  deg <- igraph::degree(g, mode = mode)
  deg <- deg[deg > 0]  # Exclude zero-degree nodes
  n <- length(deg)
  total <- sum(deg)
  p <- deg / total
  -sum(p * log2(p))
}


# ---- Inlined from /Mob/percolation_threshold.R ----
percolation_threshold <- function(g, mode = c("random", "targeted"), trials = 10) {
  mode <- match.arg(mode)
  n <- vcount(g)
  
  find_threshold <- function() {
    g_temp <- g
    for (f in seq(0.01, 1, by = 0.01)) {
      if (mode == "random") {
        remove_n <- ceiling(f * n)
        remove_v <- sample(V(g_temp), min(remove_n, vcount(g_temp)))
      } else {
        # Targeted removal: highest-degree node
        degs <- degree(g_temp)
        remove_v <- which.max(degs)
      }
      g_temp <- delete_vertices(g_temp, remove_v)
      
      if (vcount(g_temp) == 0) return(f)
      
      comps <- components(g_temp)
      lcc_frac <- max(comps$csize) / n
      
      if (lcc_frac < 0.01) return(f)  # LCC below 1% means fragmented
    }
    return(1)
  }
  
  thresholds <- replicate(trials, find_threshold())
  mean(thresholds)
}


# ---- Inlined from /Mob/graphex.R ----
#' =============================================================================
#' GRAPHEX ANALYSIS FOR DIRECTED WEIGHTED SPARSE NETWORKS
#' =============================================================================
#' 
#' This script implements:
#' 1. Graphex (Graphon + Exchangeable) matrix estimation for directed weighted graphs
#' 2. Distance computation between two graphex representations
#' 3. Pairwise comparison of multiple networks
#' 
#' Optimized for large sparse networks using sparse matrix operations.
#' 
#' Author: AI Drive Assistant
#' Date: 2025
#' =============================================================================

# -----------------------------------------------------------------------------
# 1. REQUIRED PACKAGES
# -----------------------------------------------------------------------------

required_packages <- c("Matrix", "igraph", "RSpectra", "irlba", "methods")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------------------------------------------------------
# 2. UTILITY FUNCTIONS FOR SPARSE MATRICES
# -----------------------------------------------------------------------------

#' Convert input to sparse adjacency matrix
#' 
#' Handles igraph objects, dense matrices, and sparse matrices.
#' 
#' @param x Input graph (igraph object, dense matrix, or sparse matrix)
#' @param weight_attr Name of weight attribute for igraph objects
#' @return Sparse adjacency matrix (dgCMatrix)
to_sparse_adjacency <- function(x, weight_attr = "weight") {
  
  if (inherits(x, "igraph")) {
    # Convert igraph to sparse matrix
    if (weight_attr %in% igraph::edge_attr_names(x)) {
      A <- igraph::as_adjacency_matrix(x, attr = weight_attr, sparse = TRUE)
    } else {
      A <- igraph::as_adjacency_matrix(x, sparse = TRUE)
    }
    return(as(A, "dgCMatrix"))
    
  } else if (inherits(x, "dgCMatrix") || inherits(x, "dsCMatrix")) {
    return(as(x, "dgCMatrix"))
    
  } else if (inherits(x, "sparseMatrix")) {
    return(as(x, "dgCMatrix"))
    
  } else if (is.matrix(x)) {
    # Convert dense to sparse
    return(as(x, "dgCMatrix"))
    
  } else {
    stop("Input must be an igraph object, matrix, or sparse matrix!")
  }
}


#' Compute row sums for sparse matrix efficiently
#' 
#' @param A Sparse matrix
#' @return Vector of row sums
sparse_rowSums <- function(A) {
  Matrix::rowSums(A)
}


#' Compute column sums for sparse matrix efficiently
#' 
#' @param A Sparse matrix
#' @return Vector of column sums
sparse_colSums <- function(A) {
  Matrix::colSums(A)
}


#' Count non-zero elements per row
#' 
#' @param A Sparse matrix
#' @return Vector of non-zero counts per row
sparse_row_nnz <- function(A) {
  diff(A@p)
}


#' Pad sparse matrix to target size with zeros
#' 
#' @param A Sparse matrix
#' @param target_size Target dimension (square matrix)
#' @return Padded sparse matrix
pad_sparse_matrix <- function(A, target_size) {
  current_size <- nrow(A)
  
  if (current_size >= target_size) {
    return(A[1:target_size, 1:target_size, drop = FALSE])
  }
  
  # Create larger sparse matrix
  A_new <- Matrix::sparseMatrix(
    i = A@i + 1,
    p = A@p,
    x = A@x,
    dims = c(target_size, target_size),
    dimnames = NULL
  )
  
  return(A_new)
}


# -----------------------------------------------------------------------------
# 3. GRAPHEX ESTIMATION - MAIN FUNCTION
# -----------------------------------------------------------------------------

#' Estimate Graphex matrix from directed weighted sparse graph
#' 
#' Graphex (Graphon + Exchangeable) extends graphon theory to sparse graphs
#' by incorporating degree heterogeneity and sparsity patterns.
#' 
#' The graphex representation captures:
#' - Block structure (community patterns)
#' - Degree heterogeneity (hub nodes)
#' - Sparsity patterns
#' - Directional asymmetry
#' 
#' @param A Input graph: igraph object, dense matrix, or sparse matrix
#' @param K Number of blocks/communities. If NULL, estimated automatically.
#' @param method Estimation method: "block", "spectral", or "degree_corrected"
#' @param normalize Normalization: "none", "row", "symmetric", or "laplacian"
#' @param weight_attr Weight attribute name for igraph objects
#' @return A graphex object containing the estimated representation
#' 
#' @examples
#' # Create a sparse directed weighted network
#' A <- Matrix::rsparsematrix(1000, 1000, density = 0.01)
#' A <- abs(A)  # Make weights positive
#' graphex <- estimate_graphex(A, K = 5, method = "degree_corrected")
#' 
estimate_graphex <- function(A, 
                             K = NULL, 
                             method = "degree_corrected",
                             normalize = "symmetric",
                             weight_attr = "weight") {
  
  # Convert to sparse matrix
  A <- to_sparse_adjacency(A, weight_attr)
  n <- nrow(A)
  
  # Validate input
  if (nrow(A) != ncol(A)) {
    stop("Adjacency matrix must be square!")
  }
  
  # Compute network statistics for sparsity handling
  nnz <- Matrix::nnzero(A)
  density <- nnz / (n * n)
  
  message(sprintf("Network: n=%d, edges=%d, density=%.6f", n, nnz, density))
  
  # Estimate K if not provided
  if (is.null(K)) {
    K <- estimate_K_sparse(A)
    message(sprintf("Automatically selected K=%d", K))
  }
  
  # Compute degree statistics (for degree-corrected model)
  out_strength <- sparse_rowSums(A)
  in_strength <- sparse_colSums(A)
  out_degree <- sparse_rowSums(A > 0)
  in_degree <- sparse_colSums(A > 0)
  
  # Normalize matrix based on method
  A_norm <- normalize_sparse_adjacency(A, method = normalize)
  
  # Estimate graphex based on selected method
  if (method == "block") {
    result <- estimate_graphex_block(A_norm, A, K, out_strength, in_strength)
  } else if (method == "spectral") {
    result <- estimate_graphex_spectral(A_norm, A, K)
  } else if (method == "degree_corrected") {
    result <- estimate_graphex_degree_corrected(A_norm, A, K, 
                                                out_strength, in_strength,
                                                out_degree, in_degree)
  } else {
    stop("Unknown method! Choose: 'block', 'spectral', 'degree_corrected'")
  }
  
  # Construct graphex object
  graphex <- list(
    # Core graphex components
    W = result$W,                        # Block interaction matrix (K x K)
    theta_out = result$theta_out,        # Out-degree parameters (K)
    theta_in = result$theta_in,          # In-degree parameters (K)
    rho = density,                        # Global sparsity parameter
    
    # Node assignments
    node_labels = result$labels,          # Node to block assignments
    block_sizes = as.vector(table(result$labels)),
    
    # Degree heterogeneity parameters
    degree_params = result$degree_params,
    
    # Metadata
    K = K,
    n = n,
    nnz = nnz,
    method = method,
    normalize = normalize
  )
  
  class(graphex) <- "directed_graphex"
  
  return(graphex)
}


#' Normalize sparse adjacency matrix
#' 
#' @param A Sparse adjacency matrix
#' @param method Normalization method
#' @return Normalized sparse matrix
normalize_sparse_adjacency <- function(A, method = "symmetric") {
  
  if (method == "none") {
    return(A)
  }
  
  n <- nrow(A)
  
  if (method == "row") {
    # Row normalization (stochastic matrix)
    row_sums <- sparse_rowSums(A)
    row_sums[row_sums == 0] <- 1
    D_inv <- Matrix::Diagonal(n, 1 / row_sums)
    return(D_inv %*% A)
    
  } else if (method == "symmetric") {
    # Symmetric normalization: D_out^{-1/2} A D_in^{-1/2}
    out_sums <- sparse_rowSums(A)
    in_sums <- sparse_colSums(A)
    out_sums[out_sums == 0] <- 1
    in_sums[in_sums == 0] <- 1
    D_out_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(out_sums))
    D_in_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(in_sums))
    return(D_out_inv_sqrt %*% A %*% D_in_inv_sqrt)
    
  } else if (method == "laplacian") {
    # Normalized Laplacian-style
    out_sums <- sparse_rowSums(A)
    out_sums[out_sums == 0] <- 1
    D_inv <- Matrix::Diagonal(n, 1 / out_sums)
    I <- Matrix::Diagonal(n)
    return(I - D_inv %*% A)
    
  } else {
    warning("Unknown normalization method, returning original matrix")
    return(A)
  }
}


# -----------------------------------------------------------------------------
# 4. GRAPHEX ESTIMATION METHODS
# -----------------------------------------------------------------------------

#' Estimate K (number of blocks) for sparse matrix
#' 
#' Uses eigenvalue gap heuristic on the normalized adjacency
#' 
#' @param A Sparse adjacency matrix
#' @param max_K Maximum K to consider
#' @return Estimated number of blocks
estimate_K_sparse <- function(A, max_K = 20) {
  n <- nrow(A)
  max_K <- min(max_K, floor(sqrt(n)), n - 1)
  
  # Use truncated SVD for large sparse matrices
  tryCatch({
    if (n > 100) {
      svd_result <- irlba::irlba(A, nv = max_K)
    } else {
      svd_result <- svd(as.matrix(A))
    }
    
    d <- svd_result$d[1:min(max_K, length(svd_result$d))]
    
    if (length(d) < 3) {
      return(2)
    }
    
    # Find elbow using ratio criterion
    ratios <- d[-length(d)] / (d[-1] + 1e-10)
    K <- which.max(ratios) + 1
    K <- max(2, min(K, max_K))
    
    return(K)
    
  }, error = function(e) {
    message("SVD failed, using default K=5")
    return(5)
  })
}


#' Block model estimation for graphex
#' 
#' @param A_norm Normalized sparse adjacency
#' @param A Original sparse adjacency
#' @param K Number of blocks
#' @param out_strength Out-strength vector
#' @param in_strength In-strength vector
#' @return List with W, theta_out, theta_in, labels, degree_params
estimate_graphex_block <- function(A_norm, A, K, out_strength, in_strength) {
  n <- nrow(A)
  
  # Feature matrix for clustering
  features <- cbind(
    scale_sparse_safe(out_strength),
    scale_sparse_safe(in_strength)
  )
  
  # K-means clustering
  set.seed(42)
  km <- kmeans(features, centers = K, nstart = 25, iter.max = 100)
  labels <- km$cluster
  
  # Compute block interaction matrix W
  W <- compute_block_matrix_sparse(A, labels, K)
  
  # Compute degree parameters per block
  theta_out <- compute_block_means(out_strength, labels, K)
  theta_in <- compute_block_means(in_strength, labels, K)
  
  return(list(
    W = W,
    theta_out = theta_out,
    theta_in = theta_in,
    labels = labels,
    degree_params = list(
      out_strength = out_strength,
      in_strength = in_strength
    )
  ))
}


#' Spectral estimation for graphex
#' 
#' @param A_norm Normalized sparse adjacency
#' @param A Original sparse adjacency
#' @param K Number of blocks
#' @return List with W, theta_out, theta_in, labels, degree_params
estimate_graphex_spectral <- function(A_norm, A, K) {
  n <- nrow(A)
  
  # Truncated SVD for spectral embedding
  tryCatch({
    if (n > 100) {
      svd_result <- irlba::irlba(A_norm, nv = K, nu = K)
    } else {
      svd_result <- svd(as.matrix(A_norm))
    }
    
    k_use <- min(K, length(svd_result$d))
    
    U <- svd_result$u[, 1:k_use, drop = FALSE]
    V <- svd_result$v[, 1:k_use, drop = FALSE]
    
    # Spectral embedding: concatenate left and right singular vectors
    embedding <- cbind(U, V)
    embedding[is.na(embedding)] <- 0
    embedding[!is.finite(embedding)] <- 0
    
  }, error = function(e) {
    message("SVD failed, using degree-based features")
    out_str <- sparse_rowSums(A)
    in_str <- sparse_colSums(A)
    embedding <- cbind(scale_sparse_safe(out_str), scale_sparse_safe(in_str))
  })
  
  # K-means on embedding
  set.seed(42)
  km <- kmeans(embedding, centers = K, nstart = 25, iter.max = 100)
  labels <- km$cluster
  
  # Compute block matrices and parameters
  W <- compute_block_matrix_sparse(A, labels, K)
  
  out_strength <- sparse_rowSums(A)
  in_strength <- sparse_colSums(A)
  
  theta_out <- compute_block_means(out_strength, labels, K)
  theta_in <- compute_block_means(in_strength, labels, K)
  
  return(list(
    W = W,
    theta_out = theta_out,
    theta_in = theta_in,
    labels = labels,
    degree_params = list(
      out_strength = out_strength,
      in_strength = in_strength
    )
  ))
}


#' Degree-corrected estimation for graphex
#' 
#' Implements the degree-corrected stochastic block model variant
#' which better captures hub structures in real networks.
#' 
#' @param A_norm Normalized sparse adjacency
#' @param A Original sparse adjacency
#' @param K Number of blocks
#' @param out_strength Out-strength vector
#' @param in_strength In-strength vector
#' @param out_degree Out-degree vector
#' @param in_degree In-degree vector
#' @return List with W, theta_out, theta_in, labels, degree_params
estimate_graphex_degree_corrected <- function(A_norm, A, K, 
                                              out_strength, in_strength,
                                              out_degree, in_degree) {
  n <- nrow(A)
  
  # Regularized spectral clustering for degree-corrected model
  # Use regularized Laplacian: L_tau = D_tau^{-1/2} A D_tau^{-1/2}
  # where D_tau = D + tau * I for regularization
  
  tau <- mean(c(out_degree, in_degree))  # Regularization parameter
  
  out_reg <- out_degree + tau
  in_reg <- in_degree + tau
  
  D_out_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(out_reg))
  D_in_inv_sqrt <- Matrix::Diagonal(n, 1 / sqrt(in_reg))
  
  A_reg <- D_out_inv_sqrt %*% A %*% D_in_inv_sqrt
  
  # Truncated SVD
  tryCatch({
    if (n > 100) {
      svd_result <- irlba::irlba(A_reg, nv = K, nu = K)
    } else {
      svd_result <- svd(as.matrix(A_reg))
    }
    
    k_use <- min(K, length(svd_result$d))
    
    U <- svd_result$u[, 1:k_use, drop = FALSE]
    V <- svd_result$v[, 1:k_use, drop = FALSE]
    
    # Row normalize embedding for spherical k-means effect
    U_norm <- U / (sqrt(rowSums(U^2)) + 1e-10)
    V_norm <- V / (sqrt(rowSums(V^2)) + 1e-10)
    
    embedding <- cbind(U_norm, V_norm)
    embedding[is.na(embedding)] <- 0
    embedding[!is.finite(embedding)] <- 0
    
  }, error = function(e) {
    message("SVD failed, using degree-based features")
    embedding <- cbind(
      scale_sparse_safe(out_strength),
      scale_sparse_safe(in_strength),
      scale_sparse_safe(out_degree),
      scale_sparse_safe(in_degree)
    )
  })
  
  # K-means clustering
  set.seed(42)
  km <- kmeans(embedding, centers = K, nstart = 25, iter.max = 100)
  labels <- km$cluster
  
  # Compute block matrix and degree correction parameters
  W <- compute_block_matrix_sparse(A, labels, K)
  
  # Degree correction parameters: node-level propensities
  # phi_i = d_i / sum_{j in same block} d_j
  phi_out <- numeric(n)
  phi_in <- numeric(n)
  
  for (k in 1:K) {
    idx <- which(labels == k)
    block_out_sum <- sum(out_strength[idx])
    block_in_sum <- sum(in_strength[idx])
    
    if (block_out_sum > 0) {
      phi_out[idx] <- out_strength[idx] / block_out_sum
    } else {
      phi_out[idx] <- 1 / length(idx)
    }
    
    if (block_in_sum > 0) {
      phi_in[idx] <- in_strength[idx] / block_in_sum
    } else {
      phi_in[idx] <- 1 / length(idx)
    }
  }
  
  theta_out <- compute_block_means(out_strength, labels, K)
  theta_in <- compute_block_means(in_strength, labels, K)
  
  return(list(
    W = W,
    theta_out = theta_out,
    theta_in = theta_in,
    labels = labels,
    degree_params = list(
      out_strength = out_strength,
      in_strength = in_strength,
      out_degree = out_degree,
      in_degree = in_degree,
      phi_out = phi_out,
      phi_in = phi_in
    )
  ))
}


#' Compute block interaction matrix from sparse adjacency
#' 
#' @param A Sparse adjacency matrix
#' @param labels Node to block assignments
#' @param K Number of blocks
#' @return K x K block interaction matrix
compute_block_matrix_sparse <- function(A, labels, K) {
  W <- matrix(0, K, K)
  
  for (i in 1:K) {
    idx_i <- which(labels == i)
    for (j in 1:K) {
      idx_j <- which(labels == j)
      
      if (length(idx_i) > 0 && length(idx_j) > 0) {
        block <- A[idx_i, idx_j, drop = FALSE]
        W[i, j] <- Matrix::mean(block)
      }
    }
  }
  
  return(W)
}


#' Compute mean values per block
#' 
#' @param x Numeric vector
#' @param labels Block assignments
#' @param K Number of blocks
#' @return Vector of block means
compute_block_means <- function(x, labels, K) {
  result <- numeric(K)
  for (k in 1:K) {
    idx <- which(labels == k)
    if (length(idx) > 0) {
      result[k] <- mean(x[idx])
    }
  }
  return(result)
}


#' Safe scaling for vectors with potential zero variance
#' 
#' @param x Numeric vector
#' @return Scaled vector
scale_sparse_safe <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (is.na(s) || s == 0) {
    return(x - mean(x))
  }
  return((x - mean(x)) / s)
}


# -----------------------------------------------------------------------------
# 5. GRAPHEX DISTANCE FUNCTIONS
# -----------------------------------------------------------------------------

#' Compute distance between two graphex representations
#' 
#' Implements multiple distance metrics suitable for comparing
#' sparse directed weighted networks.
#' 
#' @param graphex1 First graphex object or sparse adjacency matrix
#' @param graphex2 Second graphex object or sparse adjacency matrix
#' @param method Distance metric: "spectral", "frobenius", "block", 
#'               "wasserstein", "js", or "all"
#' @param K Number of blocks (used if inputs are matrices)
#' @return Distance value or list of distances if method="all"
#' 
#' @examples
#' G1 <- estimate_graphex(A1, K = 5)
#' G2 <- estimate_graphex(A2, K = 5)
#' d <- graphex_distance(G1, G2, method = "all")
#' 
graphex_distance <- function(graphex1, graphex2, method = "all", K = NULL) {
  
  # Extract components
  if (inherits(graphex1, "directed_graphex")) {
    W1 <- graphex1$W
    theta_out1 <- graphex1$theta_out
    theta_in1 <- graphex1$theta_in
    rho1 <- graphex1$rho
    A1 <- NULL
    n1 <- graphex1$n
  } else {
    A1 <- to_sparse_adjacency(graphex1)
    n1 <- nrow(A1)
    if (is.null(K)) K <- estimate_K_sparse(A1)
    temp <- estimate_graphex(A1, K = K)
    W1 <- temp$W
    theta_out1 <- temp$theta_out
    theta_in1 <- temp$theta_in
    rho1 <- temp$rho
  }
  
  if (inherits(graphex2, "directed_graphex")) {
    W2 <- graphex2$W
    theta_out2 <- graphex2$theta_out
    theta_in2 <- graphex2$theta_in
    rho2 <- graphex2$rho
    A2 <- NULL
    n2 <- graphex2$n
  } else {
    A2 <- to_sparse_adjacency(graphex2)
    n2 <- nrow(A2)
    if (is.null(K)) K <- estimate_K_sparse(A2)
    temp <- estimate_graphex(A2, K = K)
    W2 <- temp$W
    theta_out2 <- temp$theta_out
    theta_in2 <- temp$theta_in
    rho2 <- temp$rho
  }
  
  # Compute distances
  if (method == "all") {
    return(list(
      spectral = graphex_spectral_distance(graphex1, graphex2),
      frobenius = graphex_frobenius_distance(W1, W2),
      block = graphex_block_distance(W1, W2, theta_out1, theta_out2, 
                                     theta_in1, theta_in2),
      wasserstein = graphex_wasserstein_distance(W1, W2, theta_out1, theta_out2),
      js = graphex_js_divergence(W1, W2),
      sparsity = abs(rho1 - rho2)
    ))
  } else if (method == "spectral") {
    return(graphex_spectral_distance(graphex1, graphex2))
  } else if (method == "frobenius") {
    return(graphex_frobenius_distance(W1, W2))
  } else if (method == "block") {
    return(graphex_block_distance(W1, W2, theta_out1, theta_out2, 
                                  theta_in1, theta_in2))
  } else if (method == "wasserstein") {
    return(graphex_wasserstein_distance(W1, W2, theta_out1, theta_out2))
  } else if (method == "js") {
    return(graphex_js_divergence(W1, W2))
  } else {
    stop("Unknown method! Choose: 'spectral', 'frobenius', 'block', 'wasserstein', 'js', 'all'")
  }
}


#' Spectral distance between graphex representations
#' 
#' Computes distance based on singular value decomposition.
#' Captures global structural differences.
#' 
#' @param graphex1 First graphex or sparse matrix
#' @param graphex2 Second graphex or sparse matrix
#' @param k Number of singular values to use
#' @return Spectral distance
graphex_spectral_distance <- function(graphex1, graphex2, k = 20) {
  
  # Get adjacency matrices
  if (inherits(graphex1, "directed_graphex")) {
    W1 <- graphex1$W
    n1 <- graphex1$n
  } else {
    A1 <- to_sparse_adjacency(graphex1)
    W1 <- A1
    n1 <- nrow(A1)
  }
  
  if (inherits(graphex2, "directed_graphex")) {
    W2 <- graphex2$W
    n2 <- graphex2$n
  } else {
    A2 <- to_sparse_adjacency(graphex2)
    W2 <- A2
    n2 <- nrow(A2)
  }
  
  # Compute singular values
  k1 <- min(k, nrow(W1) - 1)
  k2 <- min(k, nrow(W2) - 1)
  
  tryCatch({
    if (nrow(W1) > 50 && inherits(W1, "sparseMatrix")) {
      svd1 <- irlba::irlba(W1, nv = k1)
    } else {
      svd1 <- svd(as.matrix(W1))
    }
    d1 <- svd1$d[1:min(k, length(svd1$d))]
  }, error = function(e) {
    d1 <- rep(0, k)
  })
  
  tryCatch({
    if (nrow(W2) > 50 && inherits(W2, "sparseMatrix")) {
      svd2 <- irlba::irlba(W2, nv = k2)
    } else {
      svd2 <- svd(as.matrix(W2))
    }
    d2 <- svd2$d[1:min(k, length(svd2$d))]
  }, error = function(e) {
    d2 <- rep(0, k)
  })
  
  # Normalize by network size
  d1 <- d1 / sqrt(n1)
  d2 <- d2 / sqrt(n2)
  
  # Align lengths
  max_len <- max(length(d1), length(d2))
  d1 <- c(d1, rep(0, max_len - length(d1)))
  d2 <- c(d2, rep(0, max_len - length(d2)))
  
  return(sqrt(sum((d1 - d2)^2)))
}


#' Frobenius distance between graphex block matrices
#' 
#' @param W1 First block matrix
#' @param W2 Second block matrix
#' @return Frobenius distance
graphex_frobenius_distance <- function(W1, W2) {
  # Align dimensions
  K <- max(nrow(W1), nrow(W2))
  W1_pad <- pad_dense_matrix(W1, K)
  W2_pad <- pad_dense_matrix(W2, K)
  
  return(norm(W1_pad - W2_pad, type = "F") / K)
}


#' Block-level distance incorporating degree parameters
#' 
#' @param W1 First block matrix
#' @param W2 Second block matrix
#' @param theta_out1 Out-degree params for first graphex
#' @param theta_out2 Out-degree params for second graphex
#' @param theta_in1 In-degree params for first graphex
#' @param theta_in2 In-degree params for second graphex
#' @return Combined block distance
graphex_block_distance <- function(W1, W2, theta_out1, theta_out2, 
                                   theta_in1, theta_in2) {
  
  # Align dimensions
  K <- max(nrow(W1), nrow(W2))
  W1_pad <- pad_dense_matrix(W1, K)
  W2_pad <- pad_dense_matrix(W2, K)
  
  # Pad theta vectors
  theta_out1 <- c(theta_out1, rep(0, K - length(theta_out1)))
  theta_out2 <- c(theta_out2, rep(0, K - length(theta_out2)))
  theta_in1 <- c(theta_in1, rep(0, K - length(theta_in1)))
  theta_in2 <- c(theta_in2, rep(0, K - length(theta_in2)))
  
  # Normalize theta vectors
  if (sum(theta_out1) > 0) theta_out1 <- theta_out1 / sum(theta_out1)
  if (sum(theta_out2) > 0) theta_out2 <- theta_out2 / sum(theta_out2)
  if (sum(theta_in1) > 0) theta_in1 <- theta_in1 / sum(theta_in1)
  if (sum(theta_in2) > 0) theta_in2 <- theta_in2 / sum(theta_in2)
  
  # Combined distance
  d_W <- norm(W1_pad - W2_pad, type = "F") / K
  d_out <- sqrt(sum((theta_out1 - theta_out2)^2))
  d_in <- sqrt(sum((theta_in1 - theta_in2)^2))
  
  return(d_W + 0.5 * (d_out + d_in))
}


#' Wasserstein distance between graphex representations
#' 
#' @param W1 First block matrix
#' @param W2 Second block matrix
#' @param theta_out1 Out-degree params for first graphex
#' @param theta_out2 Out-degree params for second graphex
#' @return Wasserstein distance approximation
graphex_wasserstein_distance <- function(W1, W2, theta_out1, theta_out2) {
  
  # Use row sums of W as distribution
  p1 <- rowSums(W1)
  p2 <- rowSums(W2)
  
  # Align lengths
  max_len <- max(length(p1), length(p2))
  p1 <- c(p1, rep(0, max_len - length(p1)))
  p2 <- c(p2, rep(0, max_len - length(p2)))
  
  # Normalize to probability distributions
  eps <- 1e-10
  p1 <- (p1 + eps) / sum(p1 + eps)
  p2 <- (p2 + eps) / sum(p2 + eps)
  
  # 1D Wasserstein (Earth Mover's Distance)
  cdf1 <- cumsum(sort(p1))
  cdf2 <- cumsum(sort(p2))
  
  return(mean(abs(cdf1 - cdf2)))
}


#' Jensen-Shannon divergence between graphex representations
#' 
#' @param W1 First block matrix
#' @param W2 Second block matrix
#' @return JS divergence
graphex_js_divergence <- function(W1, W2) {
  
  # Flatten and align
  v1 <- as.vector(W1)
  v2 <- as.vector(W2)
  
  max_len <- max(length(v1), length(v2))
  v1 <- c(v1, rep(0, max_len - length(v1)))
  v2 <- c(v2, rep(0, max_len - length(v2)))
  
  # Convert to probability distributions
  eps <- 1e-10
  v1 <- (v1 + eps) / sum(v1 + eps)
  v2 <- (v2 + eps) / sum(v2 + eps)
  
  # JS divergence
  m <- (v1 + v2) / 2
  kl1 <- sum(v1 * log(v1 / m))
  kl2 <- sum(v2 * log(v2 / m))
  
  return((kl1 + kl2) / 2)
}


#' Pad dense matrix to target size
#' 
#' @param M Dense matrix
#' @param target_size Target dimension
#' @return Padded matrix
pad_dense_matrix <- function(M, target_size) {
  current_size <- nrow(M)
  
  if (current_size >= target_size) {
    return(M[1:target_size, 1:target_size, drop = FALSE])
  }
  
  M_new <- matrix(0, target_size, target_size)
  M_new[1:current_size, 1:current_size] <- M
  
  return(M_new)
}


# -----------------------------------------------------------------------------
# 6. PAIRWISE COMPARISON OF MULTIPLE NETWORKS
# -----------------------------------------------------------------------------

#' Compare multiple networks using graphex distance
#' 
#' Computes pairwise graphex distances for a list of networks.
#' Optimized for large sparse networks.
#' 
#' @param network_list Named list of networks (igraph objects or sparse matrices)
#' @param K Number of blocks for graphex estimation (NULL for automatic)
#' @param method Graphex estimation method
#' @param distance_method Distance metric to compute
#' @param parallel Use parallel computation (requires 'parallel' package)
#' @param n_cores Number of cores for parallel computation
#' @return List containing graphex objects and distance matrices
#' 
#' @examples
#' networks <- list("2010" = A1, "2015" = A2, "2020" = A3)
#' result <- compare_graphex(networks, K = 5)
#' print(result$distances$spectral)
#' 
compare_graphex <- function(network_list, 
                            K = NULL, 
                            method = "degree_corrected",
                            distance_method = "all",
                            parallel = FALSE,
                            n_cores = 2) {
  
  n_networks <- length(network_list)
  network_names <- names(network_list)
  
  if (is.null(network_names)) {
    network_names <- paste0("Network_", 1:n_networks)
    names(network_list) <- network_names
  }
  
  message("=== Graphex Comparison ===\n")
  message(sprintf("Networks: %d", n_networks))
  
  # 1. Estimate graphex for each network
  message("\n1. Estimating graphex representations...")
  
  graphex_list <- list()
  
  for (i in 1:n_networks) {
    message(sprintf("   - %s", network_names[i]))
    graphex_list[[network_names[i]]] <- estimate_graphex(
      A = network_list[[i]],
      K = K,
      method = method
    )
  }
  
  # 2. Compute pairwise distances
  message("\n2. Computing pairwise distances...")
  
  # Initialize distance matrices
  if (distance_method == "all") {
    metrics <- c("spectral", "frobenius", "block", "wasserstein", "js", "sparsity")
  } else {
    metrics <- distance_method
  }
  
  distances <- list()
  for (m in metrics) {
    distances[[m]] <- matrix(0, n_networks, n_networks,
                             dimnames = list(network_names, network_names))
  }
  
  # Compute pairwise distances
  for (i in 1:(n_networks - 1)) {
    for (j in (i + 1):n_networks) {
      message(sprintf("   - %s vs %s", network_names[i], network_names[j]))
      
      d <- graphex_distance(
        graphex_list[[i]], 
        graphex_list[[j]], 
        method = distance_method
      )
      
      if (distance_method == "all") {
        for (m in metrics) {
          distances[[m]][i, j] <- d[[m]]
          distances[[m]][j, i] <- d[[m]]
        }
      } else {
        distances[[metrics]][i, j] <- d
        distances[[metrics]][j, i] <- d
      }
    }
  }
  
  message("\n3. Complete!")
  
  # 3. Return results
  result <- list(
    graphex_list = graphex_list,
    distances = distances,
    network_names = network_names,
    K = K,
    method = method,
    n_networks = n_networks
  )
  
  class(result) <- "graphex_comparison"
  
  return(result)
}


# -----------------------------------------------------------------------------
# 7. PRINT AND SUMMARY METHODS
# -----------------------------------------------------------------------------

#' Print method for directed_graphex objects
#' 
#' @param x A directed_graphex object
#' @param ... Additional arguments (ignored)
print.directed_graphex <- function(x, ...) {
  cat("Directed Weighted Graphex\n")
  cat("=========================\n")
  cat(sprintf("Network size (n):     %d\n", x$n))
  cat(sprintf("Number of edges:      %d\n", x$nnz))
  cat(sprintf("Density (rho):        %.6f\n", x$rho))
  cat(sprintf("Number of blocks (K): %d\n", x$K))
  cat(sprintf("Estimation method:    %s\n", x$method))
  cat(sprintf("Normalization:        %s\n", x$normalize))
  cat("\nBlock sizes:\n")
  print(x$block_sizes)
  cat("\nBlock interaction matrix (W):\n")
  print(round(x$W, 4))
  cat("\nOut-degree parameters (theta_out):\n")
  print(round(x$theta_out, 4))
  cat("\nIn-degree parameters (theta_in):\n")
  print(round(x$theta_in, 4))
  invisible(x)
}


#' Print method for graphex_comparison objects
#' 
#' @param x A graphex_comparison object
#' @param ... Additional arguments (ignored)
print.graphex_comparison <- function(x, ...) {
  cat("Graphex Comparison Results\n")
  cat("==========================\n")
  cat(sprintf("Number of networks: %d\n", x$n_networks))
  cat(sprintf("Networks: %s\n", paste(x$network_names, collapse = ", ")))
  cat(sprintf("Estimation method: %s\n", x$method))
  cat(sprintf("K: %s\n", ifelse(is.null(x$K), "automatic", x$K)))
  cat("\nAvailable distance matrices:\n")
  for (m in names(x$distances)) {
    cat(sprintf("  - %s\n", m))
  }
  cat("\nUse result$distances$<metric> to access distance matrices.\n")
  invisible(x)
}


#' Summary method for graphex_comparison
#' 
#' @param object A graphex_comparison object
#' @param ... Additional arguments (ignored)
summary.graphex_comparison <- function(object, ...) {
  cat("Graphex Comparison Summary\n")
  cat("==========================\n\n")
  
  for (m in names(object$distances)) {
    cat(sprintf("--- %s distance ---\n", m))
    D <- object$distances[[m]]
    cat("Distance matrix:\n")
    print(round(D, 4))
    
    # Compute summary statistics
    upper_tri <- D[upper.tri(D)]
    if (length(upper_tri) > 0) {
      cat(sprintf("  Min:    %.4f\n", min(upper_tri)))
      cat(sprintf("  Max:    %.4f\n", max(upper_tri)))
      cat(sprintf("  Mean:   %.4f\n", mean(upper_tri)))
      cat(sprintf("  Median: %.4f\n", median(upper_tri)))
    }
    cat("\n")
  }
  
  invisible(object)
}


# -----------------------------------------------------------------------------
# 8. EXAMPLE USAGE
# -----------------------------------------------------------------------------

#' Run example demonstrating graphex analysis
#' 
#' @return Comparison result object
run_graphex_example <- function() {
  
  cat("\n")
  cat("================================================================\n")
  cat("  EXAMPLE: Graphex Analysis for Sparse Directed Networks\n")
  cat("================================================================\n\n")
  
  # --- Generate sparse directed weighted networks ---
  
  generate_sparse_network <- function(n, density, centralization, seed) {
    set.seed(seed)
    
    # Create sparse structure
    nnz <- round(n * n * density)
    
    # Generate edges with preferential attachment-like distribution
    attractiveness <- rbeta(n, 1, 3)
    attractiveness[1:min(5, n)] <- attractiveness[1:min(5, n)] + centralization
    attractiveness <- attractiveness / sum(attractiveness)
    
    from <- sample(1:n, nnz, replace = TRUE)
    to <- sample(1:n, nnz, replace = TRUE, prob = attractiveness)
    weights <- round(rlnorm(nnz, meanlog = 1, sdlog = 1))
    
    # Remove self-loops
    valid <- from != to
    from <- from[valid]
    to <- to[valid]
    weights <- weights[valid]
    
    # Create sparse matrix
    A <- Matrix::sparseMatrix(
      i = from,
      j = to,
      x = weights,
      dims = c(n, n)
    )
    
    return(A)
  }
  
  # Generate networks
  cat("1. Generating sparse networks...\n")
  
  networks <- list(
    "2010" = generate_sparse_network(n = 500, density = 0.02, 
                                     centralization = 0.3, seed = 2010),
    "2015" = generate_sparse_network(n = 500, density = 0.025, 
                                     centralization = 0.4, seed = 2015),
    "2020" = generate_sparse_network(n = 500, density = 0.022, 
                                     centralization = 0.5, seed = 2020)
  )
  
  for (name in names(networks)) {
    A <- networks[[name]]
    cat(sprintf("   - %s: n=%d, edges=%d, density=%.4f\n", 
                name, nrow(A), Matrix::nnzero(A), 
                Matrix::nnzero(A) / (nrow(A)^2)))
  }
  
  # --- Estimate single graphex ---
  
  cat("\n2. Estimating single graphex (2020)...\n")
  
  G_2020 <- estimate_graphex(networks[["2020"]], K = 5, 
                             method = "degree_corrected")
  print(G_2020)
  
  # --- Compare two graphex ---
  
  cat("\n3. Comparing two graphex (2010 vs 2020)...\n")
  
  G_2010 <- estimate_graphex(networks[["2010"]], K = 5, 
                             method = "degree_corrected")
  
  d <- graphex_distance(G_2010, G_2020, method = "all")
  
  cat("\n   Distances:\n")
  for (m in names(d)) {
    cat(sprintf("   - %s: %.4f\n", m, d[[m]]))
  }
  
  # --- Compare all networks ---
  
  cat("\n4. Comparing all networks...\n")
  
  result <- compare_graphex(networks, K = 5, method = "degree_corrected")
  
  cat("\n   Spectral distance matrix:\n")
  print(round(result$distances$spectral, 4))
  
  cat("\n   Frobenius distance matrix:\n")
  print(round(result$distances$frobenius, 4))
  
  cat("\n================================================================\n")
  cat("  EXAMPLE COMPLETE\n")
  cat("================================================================\n")
  
  return(result)
}

# To run example:
# result <- run_graphex_example()


# ---- Inlined from /TS2VG/tsnda/R/vgprops.R ----
#-----------------------------------------------------------------------------#
#                                                                             #
#              Visibility Graph Based Forecast                                #
#                                                                             #
#  Written by: Zsolt T. Kosztyan                                              #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: August 2023                                                  #
#-----------------------------------------------------------------------------#

#' @export

vgprops <- function(ts, only_network_props=TRUE, weighted = TRUE, directed=TRUE,
                    last_n=length(ts),reverse =FALSE,
                    properties=c("SCI","SCO","SCA","DCI","DCO","DC","BC","CC",
                                 "HC","EC","AUT","HBS","PRC","AC","LC","PC",
                                 "BM","FGM","WM","DIV","LE","ECC","KNN","HUB",
                                 "SCR","Assort","DZI","DZO","DZ","BZ","CZ","EZ",
                                 "PRZ","HZ","AZ","PZ","LM","LeM","IM","SGM",
                                 "Arcs","BEZ","Dens","Diam","AVPL","ALE","GLE",
                                 "Mot","Tra","RRes","SRes","VAs","RCC"))
{
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("brainGraph", quietly = TRUE)) {
    stop(
      "Package \"brainGraph\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if ((methods::is(ts,"HVG"))|(methods::is(ts,"igraph"))){
    g<-ts
    weighted<-igraph::is.weighted(g)
    directed<-igraph::is.directed(g)
    window<-igraph::vcount(g)
  }else{
    window<-length(ts)
    g<-HVG(ts,weighted = weighted,directed = directed,reverse = reverse)
  }
  if (!igraph::is.weighted(g)) igraph::E(g)$weight<-1
  to<-window
  from<-min(max(to-last_n+1,1),to-1)

  if (igraph::is.directed(g)==TRUE){
    if (only_network_props==TRUE){
      xp<-cbind(unlist(ifelse("Assort" %in% properties,
                              igraph::assortativity.degree(g),list(NULL))),
                unlist(ifelse("DZI" %in% properties,
                              igraph::centralization.degree(g,mode="in")$centralization,list(NULL))),
                unlist(ifelse("DZO" %in% properties,
                              igraph::centralization.degree(g,mode="out")$centralization,list(NULL))),
                unlist(ifelse("DZ" %in% properties,
                              igraph::centralization.degree(g)$centralization,list(NULL))),
                unlist(ifelse("BZ" %in% properties,
                              igraph::centralization.betweenness(g)$centralization,
                              list(NULL))),
                unlist(ifelse("CZ" %in% properties,
                              igraph::centr_clo(g,mode = "all")$centralization,list(NULL))),
                unlist(ifelse("EZ" %in% properties,
                              igraph::centralization.evcent(g)$centralization,list(NULL))),
                unlist(ifelse("PRZ" %in% properties,
                              igraph::centralize(igraph::page_rank(g)$vector,
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("HZ" %in% properties,
                              igraph::centralize(igraph::harmonic_centrality(g,
                                                                             weights=abs(igraph::E(g)$weight)),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("AZ" %in% properties,
                              igraph::centralize(igraph::alpha.centrality(g,alpha=0.9),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("PZ" %in% properties,
                              igraph::centralize(igraph::power_centrality(g,exponent=0.9),
                                                 normalized=FALSE),list(NULL))),
                unlist(ifelse("LeM" %in% properties,
                              igraph::modularity(g,as.numeric(igraph::cluster_leiden(
                                igraph::as.undirected(g),
                                objective_function="modularity")$membership)),list(NULL))),
                unlist(ifelse("IM" %in% properties,
                              igraph::cluster_infomap(g)$modularity,list(NULL))),
                unlist(ifelse("SGM" %in% properties,
                              igraph::cluster_spinglass(g)$modularity,list(NULL))),
                unlist(ifelse("Arcs" %in% properties,
                              igraph::ecount(g),list(NULL))),
                unlist(ifelse("BEZ" %in% properties,
                              igraph::centralize(igraph::edge.betweenness(g,
                                                                          weights=abs(igraph::E(g)$weight)),
                                                 normalized=FALSE),list(NULL))),
                unlist(ifelse("Dens" %in% properties,
                              igraph::edge_density(g),list(NULL))),
                unlist(ifelse("Diam" %in% properties,
                              igraph::diameter(g),list(NULL))),
                unlist(ifelse("AVPL" %in% properties,
                              igraph::average.path.length(g),list(NULL))),
                unlist(ifelse("ALE" %in% properties,
                              igraph::average_local_efficiency(g),list(NULL))),
                unlist(ifelse("GLE" %in% properties,
                              igraph::global_efficiency(g),list(NULL))),
                unlist(ifelse("Mot" %in% properties,
                              igraph::count_motifs(g),list(NULL))),
                unlist(ifelse("Tra" %in% properties,
                              igraph::transitivity(g),list(NULL))),
                unlist(ifelse("RRes" %in% properties,
                              mean(percolate(g, igraph::vcount(g)/2,
                                             d = sample(igraph::V(g),
                                                        igraph::vcount(g)/2)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("SRes" %in% properties,
                              mean(percolate(g, igraph::vcount(g)/2, d = igraph::degree(g)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("VAs" %in% properties,
                              sum(abs(rowSums(as.matrix(
                                igraph::get.adjacency(g)),na.rm = TRUE)-
                                  colSums(as.matrix(igraph::get.adjacency(g)),
                                          na.rm=TRUE)))/
                                sum(abs(rowSums(as.matrix(igraph::get.adjacency(g)),
                                                na.rm = TRUE)+colSums(as.matrix(
                                                  igraph::get.adjacency(g)),na.rm=TRUE))),
                              list(NULL))),
                unlist(ifelse("RCC" %in% properties,
                              brainGraph::rich_club_coeff(g,weighted = weighted)$phi,
                              list(NULL))))
      xp[is.nan(xp)]<-0
      xp[is.na(xp)]<-0
      colnames(xp)<-c(unlist(ifelse("Assort" %in% properties,"Assort",list(NULL))),
                      unlist(ifelse("DZI" %in% properties,"DZI",list(NULL))),
                      unlist(ifelse("DZO" %in% properties,"DZO",list(NULL))),
                      unlist(ifelse("DZ" %in% properties,"DZ",list(NULL))),
                      unlist(ifelse("BZ" %in% properties,"BZ",list(NULL))),
                      unlist(ifelse("CZ" %in% properties,"CZ",list(NULL))),
                      unlist(ifelse("EZ" %in% properties,"EZ",list(NULL))),
                      unlist(ifelse("PRZ" %in% properties,"PRZ",list(NULL))),
                      unlist(ifelse("HZ" %in% properties,"HZ",list(NULL))),
                      unlist(ifelse("AZ" %in% properties,"AZ",list(NULL))),
                      unlist(ifelse("PZ" %in% properties,"PZ",list(NULL))),
                      unlist(ifelse("LeM" %in% properties,"LeM",list(NULL))),
                      unlist(ifelse("IM" %in% properties,"IM",list(NULL))),
                      unlist(ifelse("SGM" %in% properties,"SGM",list(NULL))),
                      unlist(ifelse("Arcs" %in% properties,"Arcs",list(NULL))),
                      unlist(ifelse("BEZ" %in% properties,"BEZ",list(NULL))),
                      unlist(ifelse("Dens" %in% properties,"Dens",list(NULL))),
                      unlist(ifelse("Diam" %in% properties,"Diam",list(NULL))),
                      unlist(ifelse("AVPL" %in% properties,"AVPL",list(NULL))),
                      unlist(ifelse("ALE" %in% properties,"ALE",list(NULL))),
                      unlist(ifelse("GLE" %in% properties,"GLE",list(NULL))),
                      unlist(ifelse("Mot" %in% properties,"Mot",list(NULL))),
                      unlist(ifelse("Tra" %in% properties,"Tra",list(NULL))),
                      unlist(ifelse("RRes" %in% properties,"RRes",list(NULL))),
                      unlist(ifelse("SRes" %in% properties,"SRes",list(NULL))),
                      unlist(ifelse("VAs" %in% properties,"VAs",list(NULL))),
                      unlist(ifelse("RCC" %in% properties,"RCC",list(NULL))))
    }
    else
    { #Directed, all attributes
      xp<-cbind(if ("SCI" %in% properties)
        t(igraph::graph.strength(g,mode = "in")[from:to]) else
          unlist(list(NULL)),
        if ("SCO" %in% properties)
          t(igraph::graph.strength(g,mode = "out")
            [max(from-1,1):(to-1)]) else unlist(list(NULL)),
        if ("SCA" %in% properties)
          t(igraph::graph.strength(g,mode = "all")
            [from:to]) else unlist(list(NULL)),
        if ("DCI" %in% properties)
          t(igraph::centr_degree(g,mode = "in")$res
            [from:to]) else unlist(list(NULL)),
        if ("DCO" %in% properties)
          t(igraph::centr_degree(g,mode = "out")$res
            [from:to]) else unlist(list(NULL)),
        if ("DC" %in% properties)
          t(igraph::centr_degree(g)$res
            [from:to]) else unlist(list(NULL)),
        if ("BC" %in% properties)
          t(igraph::centr_betw(g)$res
            [max(from-1,1):(to-1)]) else unlist(list(NULL)),
        if ("CC" %in% properties)
          t(igraph::closeness(g,mode="all",
                              weights=abs(igraph::E(g)$weight))
            [max(from-1,1):(to-1)]) else unlist(list(NULL)),
        if ("HC" %in% properties)
          t(igraph::harmonic_centrality(g,
                                        weights=abs(igraph::E(g)$weight))
            [max(from-1,1):(to-1)]) else unlist(list(NULL)),
        if ("EC" %in% properties)
          t(igraph::centr_eigen(g)$vector
            [from:to]) else unlist(list(NULL)),
        if ("AUT" %in% properties)
          t(igraph::authority.score(g)$vector
            [from:to]) else unlist(list(NULL)),
        if ("HBS" %in% properties)
          t(igraph::hub.score(g)$vector
            [from:to]) else unlist(list(NULL)),
        if ("PRC" %in% properties)
          t(igraph::page_rank(g)$vector
            [from:to]) else unlist(list(NULL)),
        if ("AC" %in% properties)
          t(igraph::alpha.centrality(g,alpha=0.9)
            [from:to]) else unlist(list(NULL)),
        if ("LC" %in% properties)
          t(brainGraph::centr_lev(g)
            [max(from-1,1):(to-1)]) else unlist(list(NULL)),
        if ("PC" %in% properties)
          t(igraph::power_centrality(g,exponent=0.9)
            [from:to]) else unlist(list(NULL)),
        if ("BM" %in% properties)
          t(igraph::cluster_edge_betweenness(g,
                                             weights=abs(igraph::E(g)$weight))$modularity
            [from:to]) else unlist(list(NULL)),
        if ("WM" %in% properties)
          t(igraph::cluster_walktrap(g)$modularity
            [from:to]) else unlist(list(NULL)),
        if ("LE" %in% properties)
          t(igraph::local_efficiency(g)
            [from:to]) else unlist(list(NULL)),
        if ("ECC" %in% properties)
          t(igraph::eccentricity(g)
            [from:to]) else unlist(list(NULL)),
        if ("KNN" %in% properties)
          t(igraph::graph.knn(g)$knn
            [from:to]) else unlist(list(NULL)),
        if ("HUB" %in% properties)
          t(brainGraph::hubness(g)
            [from:to]) else unlist(list(NULL)),
        if ("SCR" %in% properties)
          t(brainGraph::s_core(g)
            [from:to]) else unlist(list(NULL)),
        unlist(ifelse("Assort" %in% properties,
                      igraph::assortativity.degree(g),list(NULL))),
        unlist(ifelse("DZI" %in% properties,
                      igraph::centralization.degree(g,mode="in")$centralization,list(NULL))),
        unlist(ifelse("DZO" %in% properties,
                      igraph::centralization.degree(g,mode="out")$centralization,list(NULL))),
        unlist(ifelse("DZ" %in% properties,
                      igraph::centralization.degree(g)$centralization,list(NULL))),
        unlist(ifelse("BZ" %in% properties,
                      igraph::centralization.betweenness(g)$centralization,
                      list(NULL))),
        unlist(ifelse("CZ" %in% properties,
                      igraph::centr_clo(g,mode = "all")$centralization,list(NULL))),
        unlist(ifelse("EZ" %in% properties,
                      igraph::centralization.evcent(g)$centralization,list(NULL))),
        unlist(ifelse("PRZ" %in% properties,
                      igraph::centralize(igraph::page_rank(g)$vector,
                                         normalized = FALSE),list(NULL))),
        unlist(ifelse("HZ" %in% properties,
                      igraph::centralize(igraph::harmonic_centrality(g,
                                                                     weights=abs(igraph::E(g)$weight)),
                                         normalized = FALSE),list(NULL))),
        unlist(ifelse("AZ" %in% properties,
                      igraph::centralize(igraph::alpha.centrality(g,alpha=0.9),
                                         normalized = FALSE),list(NULL))),
        unlist(ifelse("PZ" %in% properties,
                      igraph::centralize(igraph::power_centrality(g,exponent=0.9),
                                         normalized=FALSE),list(NULL))),
        unlist(ifelse("LeM" %in% properties,
                      igraph::modularity(g,as.numeric(igraph::cluster_leiden(
                        igraph::as.undirected(g),
                        objective_function="modularity")$membership)),list(NULL))),
        unlist(ifelse("IM" %in% properties,
                      igraph::cluster_infomap(g)$modularity,list(NULL))),
        unlist(ifelse("SGM" %in% properties,
                      igraph::cluster_spinglass(g)$modularity,list(NULL))),
        unlist(ifelse("Arcs" %in% properties,
                      igraph::ecount(g),list(NULL))),
        unlist(ifelse("BEZ" %in% properties,
                      igraph::centralize(igraph::edge.betweenness(g,
                                                                  weights=abs(igraph::E(g)$weight)),
                                         normalized=FALSE),list(NULL))),
        unlist(ifelse("Dens" %in% properties,
                      igraph::edge_density(g),list(NULL))),
        unlist(ifelse("Diam" %in% properties,
                      igraph::diameter(g),list(NULL))),
        unlist(ifelse("AVPL" %in% properties,
                      igraph::average.path.length(g),list(NULL))),
        unlist(ifelse("ALE" %in% properties,
                      igraph::average_local_efficiency(g),list(NULL))),
        unlist(ifelse("GLE" %in% properties,
                      igraph::global_efficiency(g),list(NULL))),
        unlist(ifelse("Mot" %in% properties,
                      igraph::count_motifs(g),list(NULL))),
        unlist(ifelse("Tra" %in% properties,
                      igraph::transitivity(g),list(NULL))),
        unlist(ifelse("RRes" %in% properties,
                      mean(percolate(g, igraph::vcount(g)/2,
                                     d = sample(igraph::V(g),
                                                igraph::vcount(g)/2)),
                           na.rm = TRUE),list(NULL))),
        unlist(ifelse("SRes" %in% properties,
                      mean(percolate(g, igraph::vcount(g)/2,
                                     d = igraph::degree(g)),
                           na.rm = TRUE),list(NULL))),
        unlist(ifelse("VAs" %in% properties,
                      sum(abs(rowSums(as.matrix(
                        igraph::get.adjacency(g)),na.rm = TRUE)-
                          colSums(as.matrix(igraph::get.adjacency(g)),
                                  na.rm=TRUE)))/
                        sum(abs(rowSums(as.matrix(
                          igraph::get.adjacency(g)),
                          na.rm = TRUE)+colSums(as.matrix(
                            igraph::get.adjacency(g)),na.rm=TRUE))),
                      list(NULL))),
        unlist(ifelse("RCC" %in% properties,
                      brainGraph::rich_club_coeff(g,
                                                  weighted = weighted)$phi,list(NULL))))

      xp[is.nan(xp)]<-0
      xp[is.na(xp)]<-0
      colnames(xp)<-c(if ("SCI" %in% properties)
        paste("SCI",from:to,sep = "_") else unlist(list(NULL)),
        if ("SCO" %in% properties)
          paste("SCO",max(from-1,1):(to-1),sep = "_") else
            unlist(list(NULL)),
        if ("SCA" %in% properties)
          paste("SCA",from:to,sep = "_") else unlist(list(NULL)),
        if ("DCI" %in% properties)
          paste("DCI",from:to,sep = "_") else unlist(list(NULL)),
        if ("DCO" %in% properties)
          paste("DCO",from:to,sep = "_") else unlist(list(NULL)),
        if ("DC" %in% properties)
          paste("DC",from:to,sep = "_") else unlist(list(NULL)),
        if ("BC" %in% properties)
          paste("BC",max(from-1,1):(to-1),sep = "_") else
            unlist(list(NULL)),
        if ("CC" %in% properties)
          paste("CC",max(from-1,1):(to-1),sep = "_") else
            unlist(list(NULL)),
        if ("HC" %in% properties)
          paste("HC",max(from-1,1):(to-1),sep = "_") else
            unlist(list(NULL)),
        if ("EC" %in% properties)
          paste("EC",from:to,sep = "_") else unlist(list(NULL)),
        if ("AUT" %in% properties)
          paste("AUT",from:to,sep = "_") else unlist(list(NULL)),
        if ("HBS" %in% properties)
          paste("HBS",from:to,sep = "_") else unlist(list(NULL)),
        if ("PRC" %in% properties)
          paste("PRC",from:to,sep = "_") else unlist(list(NULL)),
        if ("AC" %in% properties)
          paste("AC",from:to,sep = "_") else unlist(list(NULL)),
        if ("LC" %in% properties)
          paste("LC",max(from-1,1):(to-1),sep = "_") else
            unlist(list(NULL)),
        if ("PC" %in% properties)
          paste("PC",from:to,sep = "_") else unlist(list(NULL)),
        if ("BM" %in% properties)
          paste("BM",from:to,sep = "_") else unlist(list(NULL)),
        if ("WM" %in% properties)
          paste("WM",from:to,sep = "_") else unlist(list(NULL)),
        if ("LE" %in% properties)
          paste("LE",from:to,sep = "_") else unlist(list(NULL)),
        if ("ECC" %in% properties)
          paste("ECC",from:to,sep = "_") else unlist(list(NULL)),
        if ("KNN" %in% properties)
          paste("KNN",from:to,sep = "_") else unlist(list(NULL)),
        if ("HUB" %in% properties)
          paste("HUB",from:to,sep = "_") else unlist(list(NULL)),
        if ("SCR" %in% properties)
          paste("SCR",from:to,sep = "_") else unlist(list(NULL)),
        unlist(
          ifelse("Assort" %in% properties,"Assort",list(NULL))),
        unlist(ifelse("DZI" %in% properties,"DZI",list(NULL))),
        unlist(ifelse("DZO" %in% properties,"DZO",list(NULL))),
        unlist(ifelse("DZ" %in% properties,"DZ",list(NULL))),
        unlist(ifelse("BZ" %in% properties,"BZ",list(NULL))),
        unlist(ifelse("CZ" %in% properties,"CZ",list(NULL))),
        unlist(ifelse("EZ" %in% properties,"EZ",list(NULL))),
        unlist(ifelse("PRZ" %in% properties,"PRZ",list(NULL))),
        unlist(ifelse("HZ" %in% properties,"HZ",list(NULL))),
        unlist(ifelse("AZ" %in% properties,"AZ",list(NULL))),
        unlist(ifelse("PZ" %in% properties,"PZ",list(NULL))),
        unlist(ifelse("LeM" %in% properties,"LeM",list(NULL))),
        unlist(ifelse("IM" %in% properties,"IM",list(NULL))),
        unlist(ifelse("SGM" %in% properties,"SGM",list(NULL))),
        unlist(ifelse("Arcs" %in% properties,"Arcs",list(NULL))),
        unlist(ifelse("BEZ" %in% properties,"BEZ",list(NULL))),
        unlist(ifelse("Dens" %in% properties,"Dens",list(NULL))),
        unlist(ifelse("Diam" %in% properties,"Diam",list(NULL))),
        unlist(ifelse("AVPL" %in% properties,"AVPL",list(NULL))),
        unlist(ifelse("ALE" %in% properties,"ALE",list(NULL))),
        unlist(ifelse("GLE" %in% properties,"GLE",list(NULL))),
        unlist(ifelse("Mot" %in% properties,"Mot",list(NULL))),
        unlist(ifelse("Tra" %in% properties,"Tra",list(NULL))),
        unlist(ifelse("RRes" %in% properties,"RRes",list(NULL))),
        unlist(ifelse("SRes" %in% properties,"SRes",list(NULL))),
        unlist(ifelse("VAs" %in% properties,"VAs",list(NULL))),
        unlist(ifelse("RCC" %in% properties,"RCC",list(NULL))))
    }
  }
  else
  {
    if (only_network_props==TRUE){
      xp<-cbind(unlist(ifelse("Assort" %in% properties,
                              igraph::assortativity.degree(g),list(NULL))),
                unlist(ifelse("DZ" %in% properties,
                              igraph::centralization.degree(g)$centralization,
                              list(NULL))),
                unlist(ifelse("BZ" %in% properties,
                              igraph::centralization.betweenness(g)$centralization,
                              list(NULL))),
                unlist(ifelse("CZ" %in% properties,
                              igraph::centralization.closeness(g)$centralization,
                              list(NULL))),
                unlist(ifelse("EZ" %in% properties,
                              igraph::centralization.evcent(g)$centralization,
                              list(NULL))),
                unlist(ifelse("PRZ" %in% properties,
                              igraph::centralize(igraph::page_rank(g)$vector,
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("HZ" %in% properties,
                              igraph::centralize(igraph::harmonic_centrality(g,
                                                                             weights=abs(igraph::E(g)$weight)),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("AZ" %in% properties,
                              igraph::centralize(igraph::alpha.centrality(g, alpha = 0.9),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("PZ" %in% properties,
                              igraph::centralize(igraph::power_centrality(g,exponent=0.9),
                                                 normalized=FALSE),list(NULL))),
                unlist(ifelse("LM" %in% properties,
                              mean(igraph::cluster_louvain(g)$modularity),list(NULL))),
                unlist(ifelse("LeM" %in% properties,
                              igraph::modularity(g,as.numeric(igraph::cluster_leiden(
                                igraph::as.undirected(g),
                                objective_function="modularity")$membership)),list(NULL))),
                unlist(ifelse("IM" %in% properties,
                              igraph::cluster_infomap(g)$modularity,list(NULL))),
                unlist(ifelse("SGM" %in% properties,
                              igraph::cluster_spinglass(g)$modularity,list(NULL))),
                unlist(ifelse("Arcs" %in% properties,
                              igraph::ecount(g),list(NULL))),
                unlist(ifelse("BEZ" %in% properties,
                              igraph::centralize(igraph::edge.betweenness(g,
                                                                          weights=abs(igraph::E(g)$weight)),normalized=FALSE),list(NULL))),
                unlist(ifelse("Dens" %in% properties,
                              igraph::edge_density(g),list(NULL))),
                unlist(ifelse("Diam" %in% properties,
                              igraph::diameter(g),list(NULL))),
                unlist(ifelse("AVPL" %in% properties,
                              igraph::average.path.length(g),list(NULL))),
                unlist(ifelse("ALE" %in% properties,
                              igraph::average_local_efficiency(g),list(NULL))),
                unlist(ifelse("GLE" %in% properties,
                              igraph::global_efficiency(g),list(NULL))),
                unlist(ifelse("Mot" %in% properties,
                              igraph::count_motifs(g),list(NULL))),
                unlist(ifelse("Tra" %in% properties,
                              igraph::transitivity(g),list(NULL))),
                unlist(ifelse("RRes" %in% properties,
                              mean(percolate(g, igraph::vcount(g)/2,
                                             d = sample(igraph::V(g),
                                                        igraph::vcount(g)/2)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("SRes" %in% properties,
                              mean(percolate(g, igraph::vcount(g)/2, d = igraph::degree(g)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("RCC" %in% properties,
                              brainGraph::rich_club_coeff(g,weighted = weighted)$phi,
                              list(NULL))))
      xp[is.nan(xp)]<-0
      xp[is.na(xp)]<-0
      colnames(xp)<-c(unlist(ifelse("Assort" %in% properties,"Assort",list(NULL))),
                      unlist(ifelse("DZ" %in% properties,"DZ",list(NULL))),
                      unlist(ifelse("BZ" %in% properties,"BZ",list(NULL))),
                      unlist(ifelse("CZ" %in% properties,"CZ",list(NULL))),
                      unlist(ifelse("EZ" %in% properties,"EZ",list(NULL))),
                      unlist(ifelse("PRZ" %in% properties,"PRZ",list(NULL))),
                      unlist(ifelse("HZ" %in% properties,"HZ",list(NULL))),
                      unlist(ifelse("AZ" %in% properties,"AZ",list(NULL))),
                      unlist(ifelse("PZ" %in% properties,"PZ",list(NULL))),
                      unlist(ifelse("LM" %in% properties,"LM",list(NULL))),
                      unlist(ifelse("LeM" %in% properties,"LeM",list(NULL))),
                      unlist(ifelse("IM" %in% properties,"IM",list(NULL))),
                      unlist(ifelse("SGM" %in% properties,"SGM",list(NULL))),
                      unlist(ifelse("Arcs" %in% properties,"Arcs",list(NULL))),
                      unlist(ifelse("BEZ" %in% properties,"BEZ",list(NULL))),
                      unlist(ifelse("Dens" %in% properties,"Dens",list(NULL))),
                      unlist(ifelse("Diam" %in% properties,"Diam",list(NULL))),
                      unlist(ifelse("AVPL" %in% properties,"AVPL",list(NULL))),
                      unlist(ifelse("ALE" %in% properties,"ALE",list(NULL))),
                      unlist(ifelse("GLE" %in% properties,"GLE",list(NULL))),
                      unlist(ifelse("Mot" %in% properties,"Mot",list(NULL))),
                      unlist(ifelse("Tra" %in% properties,"Tra",list(NULL))),
                      unlist(ifelse("RRes" %in% properties,"RRes",list(NULL))),
                      unlist(ifelse("SRes" %in% properties,"SRes",list(NULL))),
                      unlist(ifelse("RCC" %in% properties,"RCC",list(NULL))))
    }
    else
    {
      xp<-cbind(if ("SCA" %in% properties) t(igraph::graph.strength(g)[from:to])
                else unlist(list(NULL)),
                if ("DC" %in% properties)
                  t(igraph::centr_degree(g)$res[from:to])
                else unlist(list(NULL)),
                if ("BC" %in% properties)
                  t(igraph::centr_betw(g)$res[max(from-1,1):(to-1)])
                else unlist(list(NULL)),
                if ("CC" %in% properties)
                  t(igraph::closeness(g,
                                      weights=abs(igraph::E(g)$weight))[max(from-1,1):(to-1)])
                else unlist(list(NULL)),
                if ("HC" %in% properties)
                  t(igraph::harmonic_centrality(g,
                                                weights=abs(igraph::E(g)$weight))[max(from-1,1):(to-1)])
                else unlist(list(NULL)),
                if ("EC" %in% properties)
                  t(igraph::centr_eigen(g)$vector[from:to])
                else unlist(list(NULL)),
                if ("PRC" %in% properties)
                  t(igraph::page_rank(g)$vector[from:to])
                else unlist(list(NULL)),
                if ("AC" %in% properties)
                  t(igraph::alpha.centrality(g,alpha = 0.9)[from:to])
                else unlist(list(NULL)),
                if ("LC" %in% properties)
                  t(brainGraph::centr_lev(g)[max(from-1,1):(to-1)])
                else unlist(list(NULL)),
                if ("PC" %in% properties)
                  t(igraph::power_centrality(g,exponent=0.9)[from:to])
                else unlist(list(NULL)),
                if ("BM" %in% properties)
                  t(igraph::cluster_edge_betweenness(g,
                                                     weights=abs(igraph::E(g)$weight))$modularity[from:to])
                else unlist(list(NULL)),
                if ("FGM" %in% properties)
                  t(igraph::cluster_fast_greedy(g)$modularity[from:to])
                else unlist(list(NULL)),
                if ("WM" %in% properties)
                  t(igraph::cluster_walktrap(g)$modularity[from:to])
                else unlist(list(NULL)),
                if ("DIV" %in% properties)
                  t(igraph::graph.diversity(g)[from:to])
                else unlist(list(NULL)),
                if ("LE" %in% properties)
                  t(igraph::local_efficiency(g)[from:to])
                else unlist(list(NULL)),
                if ("ECC" %in% properties)
                  t(igraph::eccentricity(g)[from:to])
                else unlist(list(NULL)),
                if ("KNN" %in% properties)
                  t(igraph::graph.knn(g)$knn[from:to])
                else unlist(list(NULL)),
                if ("HUB" %in% properties)
                  t(brainGraph::hubness(g)[from:to])
                else unlist(list(NULL)),
                if ("SCR" %in% properties)
                  t(brainGraph::s_core(g)[from:to])
                else unlist(list(NULL)),
                unlist(ifelse("Assort" %in% properties,
                              igraph::assortativity.degree(g),list(NULL))),
                unlist(ifelse("DZ" %in% properties,
                              igraph::centralization.degree(g)$centralization,
                              list(NULL))),
                unlist(ifelse("BZ" %in% properties,
                              igraph::centralization.betweenness(g)$centralization,
                              list(NULL))),
                unlist(ifelse("CZ" %in% properties,
                              igraph::centralization.closeness(g)$centralization,
                              list(NULL))),
                unlist(ifelse("EZ" %in% properties,
                              igraph::centralization.evcent(g)$centralization,
                              list(NULL))),
                unlist(ifelse("PRZ" %in% properties,
                              igraph::centralize(igraph::page_rank(g)$vector,
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("HZ" %in% properties,
                              igraph::centralize(igraph::harmonic_centrality(g,
                                                                             weights=abs(igraph::E(g)$weight)),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("AZ" %in% properties,
                              igraph::centralize(igraph::alpha.centrality(g, alpha = 0.9),
                                                 normalized = FALSE),list(NULL))),
                unlist(ifelse("PZ" %in% properties,
                              igraph::centralize(
                                igraph::power_centrality(g,exponent=0.9),
                                normalized=FALSE),list(NULL))),
                unlist(ifelse("LM" %in% properties,
                              mean(
                                igraph::cluster_louvain(g)$modularity),list(NULL))),
                unlist(ifelse("LeM" %in% properties,
                              igraph::modularity(g,
                                                 as.numeric(igraph::cluster_leiden(
                                                   igraph::as.undirected(g),
                                                   objective_function="modularity")$membership)),list(NULL))),
                unlist(ifelse("IM" %in% properties,
                              igraph::cluster_infomap(g)$modularity,list(NULL))),
                unlist(ifelse("SGM" %in% properties,
                              igraph::cluster_spinglass(g)$modularity,list(NULL))),
                unlist(ifelse("Arcs" %in% properties,
                              igraph::ecount(g),list(NULL))),
                unlist(ifelse("BEZ" %in% properties,
                              igraph::centralize(igraph::edge.betweenness(g,
                                                                          weights=abs(igraph::E(g)$weight)),
                                                 normalized=FALSE),list(NULL))),
                unlist(ifelse("Dens" %in% properties,
                              igraph::edge_density(g),list(NULL))),
                unlist(ifelse("Diam" %in% properties,
                              igraph::diameter(g),list(NULL))),
                unlist(ifelse("AVPL" %in% properties,
                              igraph::average.path.length(g),list(NULL))),
                unlist(ifelse("ALE" %in% properties,
                              igraph::average_local_efficiency(g),list(NULL))),
                unlist(ifelse("GLE" %in% properties,
                              igraph::global_efficiency(g),list(NULL))),
                unlist(ifelse("Mot" %in% properties,
                              igraph::count_motifs(g),list(NULL))),
                unlist(ifelse("Tra" %in% properties,
                              igraph::transitivity(g),list(NULL))),
                unlist(ifelse("RRes" %in% properties,
                              mean(percolate(g, igraph::vcount(g)/2,
                                             d = sample(igraph::V(g),
                                                        igraph::vcount(g)/2)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("SRes" %in% properties,
                              mean(percolate(g,
                                             igraph::vcount(g)/2, d = igraph::degree(g)),
                                   na.rm = TRUE),list(NULL))),
                unlist(ifelse("RCC" %in% properties,
                              brainGraph::rich_club_coeff(g,
                                                          weighted = weighted)$phi,
                              list(NULL))))
      xp[is.nan(xp)]<-0
      xp[is.na(xp)]<-0
      colnames(xp)<-c(if ("SCA" %in% properties)
        paste("SC",from:to,sep = "_") else unlist(list(NULL)),
        if ("DC" %in% properties)
          paste("DC",from:to,sep = "_") else unlist(list(NULL)),
        if ("BC" %in% properties)
          paste("BC",max(from-1,1):(to-1),sep = "_")
        else unlist(list(NULL)),
        if ("CC" %in% properties)
          paste("CC",max(from-1,1):(to-1),sep = "_")
        else unlist(list(NULL)),
        if ("HC" %in% properties)
          paste("HC",max(from-1,1):(to-1),sep = "_")
        else unlist(list(NULL)),
        if ("EC" %in% properties)
          paste("EC",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("PRC" %in% properties)
          paste("PRC",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("AC" %in% properties)
          paste("AC",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("LC" %in% properties)
          paste("LC",max(from-1,1):(to-1),sep = "_")
        else unlist(list(NULL)),
        if ("PC" %in% properties)
          paste("PC",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("BM" %in% properties)
          paste("BM",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("FGM" %in% properties)
          paste("FGM",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("WM" %in% properties)
          paste("WM",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("DIV" %in% properties)
          paste("DIV",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("LE" %in% properties)
          paste("LE",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("ECC" %in% properties)
          paste("ECC",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("KNN" %in% properties)
          paste("KNN",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("HUB" %in% properties)
          paste("HUB",from:to,sep = "_")
        else unlist(list(NULL)),
        if ("SCR" %in% properties)
          paste("SCR",from:to,sep = "_")
        else unlist(list(NULL)),
        unlist(ifelse("Assort" %in% properties,"Assort",list(NULL))),
        unlist(ifelse("DZ" %in% properties,"DZ",list(NULL))),
        unlist(ifelse("BZ" %in% properties,"BZ",list(NULL))),
        unlist(ifelse("CZ" %in% properties,"CZ",list(NULL))),
        unlist(ifelse("EZ" %in% properties,"EZ",list(NULL))),
        unlist(ifelse("PRZ" %in% properties,"PRZ",list(NULL))),
        unlist(ifelse("HZ" %in% properties,"HZ",list(NULL))),
        unlist(ifelse("AZ" %in% properties,"AZ",list(NULL))),
        unlist(ifelse("PZ" %in% properties,"PZ",list(NULL))),
        unlist(ifelse("LM" %in% properties,"LM",list(NULL))),
        unlist(ifelse("LeM" %in% properties,"LeM",list(NULL))),
        unlist(ifelse("IM" %in% properties,"IM",list(NULL))),
        unlist(ifelse("SGM" %in% properties,"SGM",list(NULL))),
        unlist(ifelse("Arcs" %in% properties,"Arcs",list(NULL))),
        unlist(ifelse("BEZ" %in% properties,"BEZ",list(NULL))),
        unlist(ifelse("Dens" %in% properties,"Dens",list(NULL))),
        unlist(ifelse("Diam" %in% properties,"Diam",list(NULL))),
        unlist(ifelse("AVPL" %in% properties,"AVPL",list(NULL))),
        unlist(ifelse("ALE" %in% properties,"ALE",list(NULL))),
        unlist(ifelse("GLE" %in% properties,"GLE",list(NULL))),
        unlist(ifelse("Mot" %in% properties,"Mot",list(NULL))),
        unlist(ifelse("Tra" %in% properties,"Tra",list(NULL))),
        unlist(ifelse("RRes" %in% properties,"RRes",list(NULL))),
        unlist(ifelse("SRes" %in% properties,"SRes",list(NULL))),
        unlist(ifelse("RCC" %in% properties,"RCC",list(NULL))))
    }
  }
  class(xp)<-unique(paste(c("vgprops",class(xp))))
  return(xp)
}


# ---- Inlined from /Mob/graphon_distance_directed.R ----
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



# =============================================================================
# 1. Build all network lists from data/
# =============================================================================

log_step("Step 1: building graph lists from data/.")

DATA_DIR <- "data"
all_data_files <- list.files(DATA_DIR, full.names = FALSE)

AGGR_FILES <- grep("^[0-9]{4}_from_to_kist_weight[.]txt$", all_data_files, value = TRUE)
PE_FILES <- grep("^[0-9]{4}_from_to_kist_weight_PE[.]txt$", all_data_files, value = TRUE)
HELY_FILES <- grep("^[0-9]{4}_from_to_kist_weight_HELY_[0-9]+[.]txt$", all_data_files, value = TRUE)
SZAKT_FILES <- grep("^[0-9]{4}_from_to_kist_weight_SZAKT_[0-9]+[.]txt$", all_data_files, value = TRUE)

file_year <- function(x) substr(basename(x), 1L, 4L)
suffix_index <- function(x, tag) {
  as.integer(sub(paste0(".*_", tag, "_([0-9]+)[.]txt$"), "\\1", basename(x)))
}

AGGR_FILES <- AGGR_FILES[order(as.integer(file_year(AGGR_FILES)))]
PE_FILES <- PE_FILES[order(as.integer(file_year(PE_FILES)))]
HELY_FILES <- HELY_FILES[order(as.integer(file_year(HELY_FILES)), suffix_index(HELY_FILES, "HELY"))]
SZAKT_FILES <- SZAKT_FILES[order(as.integer(file_year(SZAKT_FILES)), suffix_index(SZAKT_FILES, "SZAKT"))]

if (!length(AGGR_FILES)) stop("No annual aggregated edge files were found in data/.")
years <- unique(file_year(AGGR_FILES))
expected_years <- as.character(2006:2024)
if (!identical(years, expected_years)) {
  stop("Aggregated inputs must cover 2006-2024 exactly. Found: ", paste(years, collapse = ", "))
}

nodes <- utils::read.table(
  file.path(DATA_DIR, "nodes.txt"), header = TRUE, stringsAsFactors = FALSE,
  check.names = FALSE
)
required_node_cols <- c("nodeID", "nodeLabel", "nodeLat", "nodeLong")
if (!all(required_node_cols %in% names(nodes))) {
  stop("data/nodes.txt must contain: ", paste(required_node_cols, collapse = ", "))
}
nodes$nodeID <- as.character(nodes$nodeID)

read_graph_file <- function(filename) {
  edges <- utils::read.table(
    file.path(DATA_DIR, filename), header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_cols <- c("from", "to", "weight")
  if (!all(required_cols %in% names(edges))) {
    stop(filename, " must contain: ", paste(required_cols, collapse = ", "))
  }
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$weight <- as.numeric(edges$weight)
  igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
}

gAGGR <- lapply(AGGR_FILES, read_graph_file)
names(gAGGR) <- file_year(AGGR_FILES)
gPE <- lapply(PE_FILES, read_graph_file)
names(gPE) <- file_year(PE_FILES)

build_layer_list <- function(files, tag) {
  out <- lapply(years, function(y) {
    f <- files[file_year(files) == y]
    if (!length(f)) return(list())
    idx <- suffix_index(f, tag)
    ord <- order(idx)
    f <- f[ord]
    idx <- idx[ord]
    if (anyDuplicated(idx)) stop("Duplicate ", tag, " index in ", y, ".")
    z <- lapply(f, read_graph_file)
    names(z) <- paste0(tag, "_", idx)
    z
  })
  names(out) <- years
  out
}
gHELY <- build_layer_list(HELY_FILES, "HELY")
gSZAKT <- build_layer_list(SZAKT_FILES, "SZAKT")

save(gAGGR, gPE, gHELY, gSZAKT, file = "calcs/graphlists.RData")
log_step("Saved calcs/graphlists.RData.")

# =============================================================================
# 2. Conventional network and focused node indicators
# =============================================================================

NETWORK_INDICATORS <- c(
  "SZI","SZO","SZA","Assort","DZI","DZO","DZ","BZ","CZ","EZ",
  "PRZ","HZ","AZ","PZ","LeM","IM","SGM","Arcs","BEZ","Dens",
  "Diam","AVPL","ALE","GLE","Mot","Tra","RRes","SRes","VAs","RCC",
  "FDIMALL","FDIMGC","KIGC","NRI","NE","PCTGC","NLAC"
)
MAIN_NETWORK_INDICATORS <- c("SZI","SZO","Dens","PRZ","LeM","Tra","NRI","Assort")
MAIN_NODE_INDICATORS <- c("SCI","SCO","PRC","COR","RC")

clean_graph <- function(g) {
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum", "ignore"))
  if (igraph::ecount(g)) {
    w <- igraph::E(g)$weight
    bad <- which(!is.finite(w))
    if (length(bad)) g <- igraph::delete_edges(g, bad)
  }
  if (!igraph::is_weighted(g)) igraph::E(g)$weight <- 1
  g
}

# Edge weights in these networks are application counts. igraph interprets the
# "weights" argument of path-based functions as EDGE LENGTH (cost), so passing
# raw counts makes a heavily used route count as long, i.e. distant. This is
# the wrong direction: a larger application flow means a stronger, closer link.
#
# Measured on the 19 yearly networks, the two conventions are not a rescaling
# of each other: the year ordering reverses (Spearman -0.40 for diameter,
# -0.41 for average path length). Diameter under the raw convention reproduced
# the previously reported 268-442 range; under the reciprocal convention it is
# below 1.
#
# Metrics that do not use path lengths (DZ, BZ, CZ, Dens, RRes, SRes) are
# bit-identical under both conventions, because igraph centr_betw() and
# centr_clo() do not accept a weights argument at all. Only Diam, AVPL, GLE,
# ALE, HZ and BEZ are affected.
#
# The raw counts are preserved in the "flow" edge attribute.
apply_edge_length_convention <- function(g, mode = c("inverse", "raw", "hop")) {
  mode <- match.arg(mode)
  if (!igraph::is_igraph(g)) return(g)
  if (!("weight" %in% igraph::edge_attr_names(g))) igraph::E(g)$weight <- 1
  if (!("flow" %in% igraph::edge_attr_names(g))) {
    igraph::E(g)$flow <- abs(as.numeric(igraph::E(g)$weight))
  }
  f <- abs(as.numeric(igraph::E(g)$flow))
  pos <- f[f > 0]
  fallback <- if (length(pos)) max(1 / pos) * 1e3 else 1
  igraph::E(g)$weight <- switch(
    mode,
    inverse = ifelse(f > 0, 1 / f, fallback),
    raw     = f,
    hop     = rep(1, length(f))
  )
  g
}

# Convention used for all path-based network indicators. Set to "raw" only to
# reproduce the pre-correction numbers.
EDGE_LENGTH_MODE <- Sys.getenv("MOB_EDGE_LENGTH_MODE", "inverse")

giant_weak_component_impl <- function(g) {
  cmp <- igraph::components(g, mode = "weak")
  if (!length(cmp$csize)) return(g)
  igraph::induced_subgraph(g, igraph::V(g)[cmp$membership == which.max(cmp$csize)])
}
giant_weak_component <- giant_weak_component_impl

calc_network_indicators <- function(g, context = "network") {
  g <- clean_graph(g)
  GC <- giant_weak_component(g)
  n <- igraph::vcount(g)

  s_in <- igraph::strength(g, mode = "in")
  s_out <- igraph::strength(g, mode = "out")
  s_all <- igraph::strength(g, mode = "all")

  SZI <- igraph::centralize(s_in, max(s_in * n))
  SZO <- igraph::centralize(s_out, max(s_out * n))
  SZA <- igraph::centralize(s_all, max(s_all * n))

  # Author's original tsnda::vgprops() network-level logic is embedded above.
  # Path-based indicators need edge LENGTHS, not raw application counts.
  GC_len <- apply_edge_length_convention(GC, EDGE_LENGTH_MODE)
  NET <- tryCatch(vgprops(GC_len, only_network_props = TRUE), error = function(e) {
    stop(context, ": vgprops failed: ", conditionMessage(e), call. = FALSE)
  })
  NET <- as.numeric(NET)
  if (length(NET) != 27L) {
    stop(context, ": vgprops returned ", length(NET), " values; expected 27.")
  }
  names(NET) <- c(
    "Assort","DZI","DZO","DZ","BZ","CZ","EZ","PRZ","HZ","AZ","PZ",
    "LeM","IM","SGM","Arcs","BEZ","Dens","Diam","AVPL","ALE","GLE",
    "Mot","Tra","RRes","SRes","VAs","RCC"
  )

  FDIMALL <- safe_scalar(network_box_counting(g)$dimension, paste0(context, " FDIMALL"))
  FDIMGC <- safe_scalar(network_box_counting(GC)$dimension, paste0(context, " FDIMGC"))
  KIGC <- safe_scalar(kirchhoff_index(GC), paste0(context, " KIGC"))
  NRI <- safe_scalar(normalized_resilience_index(g), paste0(context, " NRI"))
  NE <- safe_scalar(network_entropy(g), paste0(context, " NE"))
  PCTGC <- safe_scalar(percolation_threshold(GC), paste0(context, " PCTGC"))
  NLAC <- round(mean(replicate(NLAC_REPS, mean(network_lacunarity(g)$lacunarity))), 3)

  out <- c(SZI = SZI, SZO = SZO, SZA = SZA, NET,
           FDIMALL = FDIMALL, FDIMGC = FDIMGC, KIGC = KIGC,
           NRI = NRI, NE = NE, PCTGC = PCTGC, NLAC = NLAC)
  out <- out[NETWORK_INDICATORS]
  if (length(out) != 37L || anyDuplicated(names(out))) {
    stop(context, ": invalid 37-indicator result.")
  }
  out
}

calc_main_node_indicators <- function(g) {
  g <- clean_graph(g)
  vals <- rbind(
    SCI = igraph::strength(g, mode = "in"),
    SCO = igraph::strength(g, mode = "out"),
    PRC = igraph::page_rank(g)$vector,
    COR = igraph::coreness(g, mode = "all"),
    RC = resilience_centrality(g)$resilience_centrality
  )
  vals[!is.finite(vals)] <- NA_real_
  vals
}

calc_series_matrix <- function(graphs, prefix) {
  out <- matrix(NA_real_, nrow = length(graphs), ncol = length(NETWORK_INDICATORS),
                dimnames = list(names(graphs), NETWORK_INDICATORS))
  for (i in seq_along(graphs)) {
    log_step("Network metrics: ", prefix, " ", names(graphs)[i], ".")
    out[i, ] <- calc_network_indicators(graphs[[i]], paste(prefix, names(graphs)[i]))
  }
  out
}

calc_layer_array <- function(graphs_by_year, n_layers, prefix) {
  out <- array(
    NA_real_, dim = c(length(NETWORK_INDICATORS), n_layers, length(years)),
    dimnames = list(NETWORK_INDICATORS, as.character(seq_len(n_layers)), years)
  )
  for (iy in seq_along(years)) {
    gy <- graphs_by_year[[iy]]
    if (!length(gy)) next
    for (j in seq_along(gy)) {
      if (j > n_layers) stop(prefix, " has more than ", n_layers, " layers in ", years[iy], ".")
      log_step("Network metrics: ", prefix, " ", years[iy], " layer ", j, ".")
      out[, j, iy] <- calc_network_indicators(
        gy[[j]], paste(prefix, years[iy], names(gy)[j])
      )
    }
  }
  out
}

set.seed(SEED)
AGGRNETMTX <- calc_series_matrix(gAGGR, "AGGR")
PENETMTX <- calc_series_matrix(gPE, "PE")
NET3DHELY <- calc_layer_array(gHELY, 3L, "HELY")
NET3DSZAKT <- calc_layer_array(gSZAKT, 13L, "SZAKT")

NOD3DAGGR <- array(
  NA_real_, dim = c(length(MAIN_NODE_INDICATORS), nrow(nodes), length(years)),
  dimnames = list(MAIN_NODE_INDICATORS, as.character(nodes$nodeID), years)
)
for (i in seq_along(gAGGR)) {
  vals <- calc_main_node_indicators(gAGGR[[i]])
  ids <- as.character(igraph::V(gAGGR[[i]])$name)
  idx <- match(ids, nodes$nodeID)
  if (anyNA(idx)) stop("Node alignment failed for ", years[i], ".")
  NOD3DAGGR[, idx, i] <- vals
}
AGGRNODMTX <- do.call(rbind, lapply(seq_along(years), function(i) {
  x <- t(NOD3DAGGR[, , i, drop = FALSE][, , 1L])
  colnames(x) <- MAIN_NODE_INDICATORS
  as.data.frame(x)
}))
save(AGGRNETMTX, PENETMTX, NET3DHELY, NET3DSZAKT, file = "calcs/netprops.RData")
save(AGGRNODMTX, NOD3DAGGR, file = "calcs/nodeprops.RData")
log_step("Saved calcs/netprops.RData and calcs/nodeprops.RData.")

# =============================================================================
# 3. Conventional-metric audit, NDR objects and Excel figure-data outputs
# =============================================================================

audit_dir <- "output/revision/metric_audit"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

coverage <- data.frame(
  indicator = colnames(AGGRNETMTX),
  finite_n = vapply(as.data.frame(AGGRNETMTX), function(x) sum(is.finite(x)), integer(1)),
  missing_n = vapply(as.data.frame(AGGRNETMTX), function(x) sum(!is.finite(x)), integer(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(coverage, file.path(audit_dir, "METRIC_AUDIT_coverage.csv"), row.names = FALSE)

rho <- suppressWarnings(stats::cor(AGGRNETMTX, use = "pairwise.complete.obs", method = "spearman"))
utils::write.csv(cbind(indicator = rownames(rho), as.data.frame(rho, check.names = FALSE)),
                 file.path(audit_dir, "METRIC_AUDIT_spearman.csv"), row.names = FALSE)

rho_cluster <- rho
rho_cluster[!is.finite(rho_cluster)] <- 0
diag(rho_cluster) <- 1
hc <- stats::hclust(stats::as.dist(1 - abs(rho_cluster)), method = "complete")
clusters <- data.frame(
  indicator = names(stats::cutree(hc, h = 0.10)),
  redundancy_cluster = as.integer(stats::cutree(hc, h = 0.10)),
  threshold = 0.90,
  stringsAsFactors = FALSE
)
clusters <- clusters[order(clusters$redundancy_cluster, clusters$indicator), ]
utils::write.csv(clusters, file.path(audit_dir, "METRIC_AUDIT_redundancy_clusters.csv"),
                 row.names = FALSE)

pca_x <- as.data.frame(AGGRNETMTX)
keep_pca <- vapply(pca_x, function(x) {
  z <- x[is.finite(x)]
  length(z) >= 3L && stats::sd(z) > 0
}, logical(1))
pca_x <- pca_x[, keep_pca, drop = FALSE]
for (nm in names(pca_x)) {
  med <- stats::median(pca_x[[nm]][is.finite(pca_x[[nm]])], na.rm = TRUE)
  pca_x[[nm]][!is.finite(pca_x[[nm]])] <- med
}
pca_fit <- stats::prcomp(pca_x, center = TRUE, scale. = TRUE)
pca_variance <- data.frame(
  component = paste0("PC", seq_along(pca_fit$sdev)),
  variance_share = pca_fit$sdev^2 / sum(pca_fit$sdev^2),
  cumulative_share = cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2)),
  stringsAsFactors = FALSE
)
utils::write.csv(pca_variance, file.path(audit_dir, "METRIC_AUDIT_pca_variance.csv"),
                 row.names = FALSE)
utils::write.csv(cbind(indicator = rownames(pca_fit$rotation), pca_fit$rotation),
                 file.path(audit_dir, "METRIC_AUDIT_pca_loadings.csv"), row.names = FALSE)

set.seed(1)
ORIG <- AGGRNETMTX
NDR_ORIG <- nda::ndr(ORIG, cor_method = 4)
ORIG_1 <- t(NET3DHELY[, 1, ])
ORIG_2 <- t(NET3DHELY[, 2, ])
ORIG_3 <- t(NET3DHELY[, 3, ])
NDR_1 <- nda::ndr(ORIG_1, cor_method = 4)
NDR_2 <- nda::ndr(ORIG_2, cor_method = 4)
NDR_3 <- nda::ndr(ORIG_3, cor_method = 4)
save(ORIG, NDR_ORIG, ORIG_1, ORIG_2, ORIG_3, NDR_1, NDR_2, NDR_3,
     file = "calcs/ndr_objects.RData")

professions <- c(
  "Agriculture","Human Studies","Social Sciences","Information Technology",
  "Law and Public Administration","Military","Business and Economics",
  "Engineering","Health and Medical Sciences","Pedagogy","Sport Sciences",
  "Natural Sciences","Arts"
)

safe_sheet <- function(x, used = character()) {
  x <- gsub("[\\\\/?*\\[\\]:]", "_", x)
  x <- substr(x, 1L, 31L)
  if (!x %in% used) return(x)
  k <- 2L
  repeat {
    y <- paste0(substr(x, 1L, 27L), "_", k)
    if (!y %in% used) return(y)
    k <- k + 1L
  }
}
add_sheet_data <- function(wb, name, x, used) {
  nm <- safe_sheet(name, used)
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, as.data.frame(x, check.names = FALSE), rowNames = TRUE)
  c(used, nm)
}

fig01 <- data.frame(Professions = professions, stringsAsFactors = FALSE, check.names = FALSE)
for (iy in seq_along(years)) {
  fig01[[years[iy]]] <- vapply(seq_len(13L), function(j) {
    g <- gSZAKT[[iy]][[j]]
    if (is.null(g) || !igraph::ecount(g)) return(0)
    sum(as.numeric(igraph::E(g)$weight), na.rm = TRUE)
  }, numeric(1))
}
openxlsx::write.xlsx(fig01, "output/FIG01.xlsx", overwrite = TRUE)

main_years <- c("2006","2015","2024")
tab_main <- as.data.frame(t(AGGRNETMTX[main_years, MAIN_NETWORK_INDICATORS, drop = FALSE]),
                          check.names = FALSE)
tab_main$Indicator <- rownames(tab_main)
tab_main <- tab_main[, c("Indicator", main_years)]
rownames(tab_main) <- NULL
for (j in main_years) tab_main[[j]] <- round(as.numeric(tab_main[[j]]), 4)
openxlsx::write.xlsx(tab_main, "output/FIG02_MAIN_NETWORK_INDICATORS.xlsx",
                     sheetName = "Selected indicators", overwrite = TRUE)

z <- as.data.frame(scale(ORIG[, MAIN_NETWORK_INDICATORS, drop = FALSE]), check.names = FALSE)
z$Year <- as.integer(rownames(ORIG))
openxlsx::write.xlsx(z, "output/FIG03_FOCUSED_NETWORK_TRAJECTORIES.xlsx",
                     sheetName = "z scores", overwrite = TRUE)

node_years <- c("2006","2024")
regions <- c(Budapest = "24", Debrecen = "33", Szeged = "140", Pecs = "118")
node_tabs <- lapply(names(regions), function(region_name) {
  a <- NOD3DAGGR[MAIN_NODE_INDICATORS, regions[[region_name]], node_years, drop = FALSE]
  m <- matrix(as.numeric(a), nrow = length(MAIN_NODE_INDICATORS),
              ncol = length(node_years), dimnames = list(MAIN_NODE_INDICATORS, node_years))
  ans <- data.frame(Region = region_name, Indicator = MAIN_NODE_INDICATORS,
                    Y2006 = round(m[, "2006"], 4), Y2024 = round(m[, "2024"], 4),
                    check.names = FALSE)
  names(ans)[3:4] <- c("2006", "2024")
  ans
})
tab_node <- do.call(rbind, node_tabs)
rownames(tab_node) <- NULL
openxlsx::write.xlsx(tab_node, "output/FIG03_MAIN_NODE_INDICATORS.xlsx",
                     sheetName = "Selected node indicators", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
used <- character()
used <- add_sheet_data(wb, "Network-Level Indicators", ORIG, used)
used <- add_sheet_data(wb, "Loadings (Net.Ind, GNDA)", NDR_ORIG$loadings, used)
used <- add_sheet_data(wb, "Scores (Net.Ind, GNDA)", NDR_ORIG$scores, used)
openxlsx::saveWorkbook(wb, "output/FIGA1_NETIND_ALL.xlsx", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
used <- character()
for (entry in list(
  list("Indicators (1st pref)", ORIG_1), list("Loadings (Ind., GNDA, 1st)", NDR_1$loadings),
  list("Scores (Ind., GNDA, 1st)", NDR_1$scores), list("Indicators (2nd pref)", ORIG_2),
  list("Loadings (Ind., 2nd)", NDR_2$loadings), list("Scores (Ind., 2nd)", NDR_2$scores),
  list("Indicators (3rd pref)", ORIG_3), list("Loadings (Ind., 3rd)", NDR_3$loadings),
  list("Scores (Ind., 3rd)", NDR_3$scores)
)) used <- add_sheet_data(wb, entry[[1]], entry[[2]], used)
openxlsx::saveWorkbook(wb, "output/FIGA2_PREFERENCES.xlsx", overwrite = TRUE)

wb <- openxlsx::createWorkbook()
used <- character()
for (i in seq_len(13L)) {
  DATAP <- t(NET3DSZAKT[, i, ])
  NDR_P <- nda::ndr(DATAP, cor_method = 4)
  used <- add_sheet_data(wb, professions[i], DATAP, used)
  used <- add_sheet_data(wb, paste("Loadings", professions[i]), NDR_P$loadings, used)
  used <- add_sheet_data(wb, paste("Scores", professions[i]), NDR_P$scores, used)
}
openxlsx::saveWorkbook(wb, "output/FIGA3_PROFESSIONS.xlsx", overwrite = TRUE)
log_step("Saved metric audit and manuscript Excel tables.")

# =============================================================================
# 4. Legacy graphon objects retained for Mob_final.Rmd compatibility
# =============================================================================

if (!SKIP_LEGACY_GRAPHONS) {
  log_step("Step 4: computing legacy graphon-distance objects.")
  set.seed(SEED)
  dAGGR <- compare_graphons(gAGGR)
  dPE <- compare_graphons(gPE)
  dHELY_EV <- lapply(gHELY, function(x) if (length(x) >= 2L) compare_graphons(x) else NULL)
  names(dHELY_EV) <- names(gHELY)

  hely_ids <- sort(unique(unlist(lapply(gHELY, names))))
  dEV_HELY <- lapply(hely_ids, function(id) {
    z <- lapply(gHELY, function(x) x[[id]])
    keep <- !vapply(z, is.null, logical(1))
    z <- z[keep]
    names(z) <- names(gHELY)[keep]
    if (length(z) >= 2L) compare_graphons(z) else NULL
  })
  names(dEV_HELY) <- hely_ids

  dSZAKT_EV <- lapply(gSZAKT, function(x) if (length(x) >= 2L) compare_graphons(x) else NULL)
  names(dSZAKT_EV) <- names(gSZAKT)
  szakt_ids <- unique(unlist(lapply(gSZAKT, names)))
  szakt_ids <- szakt_ids[order(as.integer(sub("SZAKT_", "", szakt_ids)))]
  dEV_SZAKT <- lapply(szakt_ids, function(id) {
    z <- lapply(gSZAKT, function(x) x[[id]])
    keep <- !vapply(z, is.null, logical(1))
    z <- z[keep]
    names(z) <- names(gSZAKT)[keep]
    if (length(z) >= 2L) compare_graphons(z) else NULL
  })
  names(dEV_SZAKT) <- szakt_ids
} else {
  warning("Legacy graphon objects were skipped by MOB_SKIP_LEGACY_GRAPHONS=1.")
  dAGGR <- dPE <- dHELY_EV <- dEV_HELY <- dSZAKT_EV <- dEV_SZAKT <- list()
}
save(dAGGR, dPE, dHELY_EV, dEV_HELY, dSZAKT_EV, dEV_SZAKT,
     file = "calcs/graphon_dists.RData")

# =============================================================================
# 5-8. Reviewer-aligned graphon, role, multi-resolution and flow modules
# =============================================================================

local({
# =============================================================================
# run_for_rev_final.R
# Final reviewer-aligned graphon robustness outputs for the Mob revision
# =============================================================================
options(stringsAsFactors = FALSE)

find_project_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/", mustWork = TRUE)))
  }
  if (file.exists("Mob.Rproj") || file.exists("Mob.Rmd")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  stop("Cannot locate the Mob project root. Run this script from /Mob.")
}

ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(ROOT)
SEED <- 20260811L
set.seed(SEED)
K_REF <- 13L
K_GRID <- c(2L, 3L, 4L, 5L, 6L, 8L, 10L, 13L)
METHODS <- c("block", "smooth", "spectral")
GRAPHON_METRICS <- c("frobenius", "cut", "wasserstein", "jensen_shannon", "total_variation")
PRIMARY_METRICS <- c("frobenius", "cut")
SECONDARY_METRICS <- c("wasserstein")
EXPLORATORY_METRICS <- c("jensen_shannon", "total_variation")

req <- c("calcs/graphlists.RData", "graphon_distance_directed.R")
miss <- req[!file.exists(req)]
if (length(miss)) stop("Missing required file(s): ", paste(miss, collapse = ", "))

pkgs <- c("igraph", "Matrix", "RSpectra", "ggplot2", "reshape2")
missp <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missp)) stop("Missing required R package(s): ", paste(missp, collapse = ", "))

suppressPackageStartupMessages(library(igraph))
# Inlined graphon/geography functions are already available.
load("calcs/graphlists.RData")

dir.create("output/revision", recursive = TRUE, showWarnings = FALSE)
dir.create("calcs/revision", recursive = TRUE, showWarnings = FALSE)

safe_spearman <- function(x, y) {
  k <- is.finite(x) & is.finite(y)
  if (sum(k) < 3L || length(unique(x[k])) < 2L || length(unique(y[k])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[k], y[k], method = "spearman"))
}

as_consecutive_series <- function(result, metric) {
  D <- result$distances[[metric]]
  if (is.null(D) || nrow(D) < 2L) return(NULL)
  x <- vapply(seq_len(nrow(D) - 1L), function(i) D[i, i + 1L], numeric(1))
  y <- result$network_names
  names(x) <- paste0(y[-length(y)], "-", y[-1L])
  x
}

rank_top_years <- function(x, n = 5L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(character())
  names(sort(x, decreasing = TRUE))[seq_len(min(n, length(x)))]
}

jaccard <- function(a, b) {
  u <- union(unique(a), unique(b))
  if (!length(u)) return(NA_real_)
  length(intersect(unique(a), unique(b))) / length(u)
}

cat("Final graphon revision analysis\n")
cat("Reference resolution: K =", K_REF, "\n")

# -----------------------------------------------------------------------------
# 1. Feasibility and specification
# -----------------------------------------------------------------------------
feasible_by_year <- vapply(
  gAGGR,
  function(g) n_distinct_rows(block_features(as_weight_matrix(g))),
  integer(1)
)

if (any(feasible_by_year < K_REF)) {
  bad <- names(feasible_by_year)[feasible_by_year < K_REF]
  stop("K=13 is not feasible for: ", paste(bad, collapse = ", "))
}

specification <- data.frame(
  item = c(
    "reference_resolution", "reference_estimator", "primary_distance_1",
    "primary_distance_2", "secondary_distance", "exploratory_distances",
    "weighting", "comparison_surface", "spectral_benchmark", "adaptive_K_role"
  ),
  value = c(
    as.character(K_REF), "block", "frobenius", "cut", "wasserstein",
    "jensen_shannon; total_variation", "raw annual application-count weights",
    "node-aligned 175x175 fitted kernels",
    "singular-spectrum distance on original weighted adjacency matrices",
    "diagnostic only; not interpreted as structural complexity"
  ),
  stringsAsFactors = FALSE
)
write.csv(specification, "output/revision/FINAL_graphon_specification.csv", row.names = FALSE)

feasibility <- data.frame(
  year = names(feasible_by_year),
  distinct_block_feature_profiles = as.integer(feasible_by_year),
  K_reference = K_REF,
  feasible = as.integer(feasible_by_year) >= K_REF,
  stringsAsFactors = FALSE
)
write.csv(feasibility, "output/revision/FINAL_K13_feasibility.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Estimator robustness at the SAME fixed K=13
# -----------------------------------------------------------------------------
cat("Estimating block/smooth/spectral graphons at fixed K=13 ...\n")
fixed_results <- setNames(vector("list", length(METHODS)), METHODS)
for (m in METHODS) {
  set.seed(SEED)
  fixed_results[[m]] <- compare_graphons(gAGGR, K = K_REF, method = m)
}

trajectory_rows <- list()
r <- 1L
for (m in METHODS) {
  for (metric in GRAPHON_METRICS) {
    x <- as_consecutive_series(fixed_results[[m]], metric)
    trajectory_rows[[r]] <- data.frame(
      estimator = m,
      fixed_K = K_REF,
      metric = metric,
      year_pair = names(x),
      distance = as.numeric(x),
      metric_role = ifelse(
        metric %in% PRIMARY_METRICS, "primary",
        ifelse(metric %in% SECONDARY_METRICS, "secondary", "exploratory")
      ),
      stringsAsFactors = FALSE
    )
    r <- r + 1L
  }
}
fixed_trajectories <- do.call(rbind, trajectory_rows)
write.csv(
  fixed_trajectories,
  "output/revision/FINAL_fixedK13_estimator_trajectories.csv",
  row.names = FALSE
)

agreement_rows <- list()
r <- 1L
for (metric in GRAPHON_METRICS) {
  for (i in seq_len(length(METHODS) - 1L)) {
    for (j in seq.int(i + 1L, length(METHODS))) {
      a <- as_consecutive_series(fixed_results[[METHODS[i]]], metric)
      b <- as_consecutive_series(fixed_results[[METHODS[j]]], metric)
      co <- intersect(names(a), names(b))
      agreement_rows[[r]] <- data.frame(
        fixed_K = K_REF,
        metric = metric,
        metric_role = ifelse(
          metric %in% PRIMARY_METRICS, "primary",
          ifelse(metric %in% SECONDARY_METRICS, "secondary", "exploratory")
        ),
        estimator_1 = METHODS[i],
        estimator_2 = METHODS[j],
        spearman_rho = safe_spearman(a[co], b[co]),
        top5_jaccard = jaccard(rank_top_years(a), rank_top_years(b)),
        stringsAsFactors = FALSE
      )
      r <- r + 1L
    }
  }
}
fixed_estimator_agreement <- do.call(rbind, agreement_rows)
write.csv(
  fixed_estimator_agreement,
  "output/revision/FINAL_fixedK13_estimator_agreement.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Fixed-K sensitivity, now compared with the fixed-K=13 reference
# -----------------------------------------------------------------------------
cat("Computing fixed-K sensitivity relative to K=13 ...\n")
k_results <- setNames(vector("list", length(K_GRID)), paste0("K", K_GRID))
for (k in K_GRID) {
  if (any(feasible_by_year < k)) {
    stop("Requested K=", k, " is not feasible in every year.")
  }
  set.seed(SEED)
  k_results[[paste0("K", k)]] <- compare_graphons(gAGGR, K = k, method = "block")
}

ksens_trajectory_rows <- list()
r <- 1L
for (k in K_GRID) {
  for (metric in GRAPHON_METRICS) {
    x <- as_consecutive_series(k_results[[paste0("K", k)]], metric)
    ksens_trajectory_rows[[r]] <- data.frame(
      metric = metric,
      metric_role = ifelse(
        metric %in% PRIMARY_METRICS, "primary",
        ifelse(metric %in% SECONDARY_METRICS, "secondary", "exploratory")
      ),
      fixed_K = k,
      year_pair = names(x),
      distance = as.numeric(x),
      stringsAsFactors = FALSE
    )
    r <- r + 1L
  }
}
ksens_trajectories <- do.call(rbind, ksens_trajectory_rows)
write.csv(
  ksens_trajectories,
  "output/revision/FINAL_K_sensitivity_trajectories.csv",
  row.names = FALSE
)

ksens_agreement_rows <- list()
r <- 1L
for (metric in GRAPHON_METRICS) {
  ref <- as_consecutive_series(k_results[[paste0("K", K_REF)]], metric)
  for (k in K_GRID) {
    x <- as_consecutive_series(k_results[[paste0("K", k)]], metric)
    co <- intersect(names(ref), names(x))
    ksens_agreement_rows[[r]] <- data.frame(
      metric = metric,
      metric_role = ifelse(
        metric %in% PRIMARY_METRICS, "primary",
        ifelse(metric %in% SECONDARY_METRICS, "secondary", "exploratory")
      ),
      fixed_K = k,
      reference_K = K_REF,
      spearman_rho_vs_K13 = safe_spearman(ref[co], x[co]),
      top5_jaccard_vs_K13 = jaccard(rank_top_years(ref), rank_top_years(x)),
      stringsAsFactors = FALSE
    )
    r <- r + 1L
  }
}
ksens_agreement <- do.call(rbind, ksens_agreement_rows)
write.csv(
  ksens_agreement,
  "output/revision/FINAL_K_sensitivity_vs_K13.csv",
  row.names = FALSE
)

# A compact main-text transition ranking for the reference block estimator.
ref_block <- fixed_trajectories[
  fixed_trajectories$estimator == "block" &
    fixed_trajectories$metric %in% c(PRIMARY_METRICS, SECONDARY_METRICS),
]
ref_block$rank_within_metric <- ave(
  -ref_block$distance,
  ref_block$metric,
  FUN = function(z) rank(z, ties.method = "min")
)
ref_block <- ref_block[order(ref_block$metric, ref_block$rank_within_metric), ]
write.csv(
  ref_block,
  "output/revision/FINAL_primary_transition_ranking.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Separate adjacency singular-spectrum benchmark
# -----------------------------------------------------------------------------
bench <- numeric(length(gAGGR) - 1L)
years <- names(gAGGR)
for (i in seq_along(bench)) {
  bench[i] <- singular_spectrum_distance(
    as_weight_matrix(gAGGR[[i]]),
    as_weight_matrix(gAGGR[[i + 1L]])
  )
}
singular_benchmark <- data.frame(
  year_pair = paste0(years[-length(years)], "-", years[-1L]),
  distance = bench,
  benchmark = "singular spectrum of original weighted adjacency matrices",
  stringsAsFactors = FALSE
)
write.csv(
  singular_benchmark,
  "output/revision/FINAL_singular_spectrum_benchmark.csv",
  row.names = FALSE
)

save(
  K_REF, K_GRID, specification, feasibility,
  fixed_results, fixed_trajectories, fixed_estimator_agreement,
  k_results, ksens_trajectories, ksens_agreement,
  singular_benchmark,
  file = "calcs/revision/FINAL_graphon_revision_objects.RData"
)

writeLines(capture.output(sessionInfo()), "output/revision/FINAL_sessionInfo.txt")

cat("\nFinal reviewer-aligned graphon analysis completed.\n")
cat("Main reference: fixed K=13 block estimator.\n")
cat("Primary distances: Frobenius and cut.\n")
cat("Secondary distance: Wasserstein.\n")
cat("JS/TV: exploratory/supplementary only.\n")
cat("Singular spectrum: separate adjacency benchmark.\n")
cat("No manuscript or existing input file was modified.\n")


})

local({
# =============================================================================
# run_graphon_fixedK13_validation.R
# Fixed-K=13 structural-role validation for the Mob revision
# =============================================================================
#
# Scientific purpose
#   - Validate whether differentiation among destination-oriented higher-
#     education regions persists when block resolution is held fixed at K=13.
#   - Focus on 2015, 2018, 2019, 2020, 2021, and 2024 so that the late-2010s
#     structural break can be inspected without conflating it with adaptive K.
#   - Avoid treating raw block labels as temporally stable. Cross-year summaries
#     use node-aligned fitted roles, within-year block intensity ranks, overlap
#     tables, and label-invariant Adjusted Rand Index (ARI).
#   - Define a destination-oriented reference core transparently from the data:
#     the top 20 micro-regions by mean incoming application weight over the full
#     2006-2024 AGGR series. This reference set is descriptive only and is not
#     used to fit the graphon/block model.
#
# This script does NOT modify Mob.Rmd, graphlists, graphon distances, or any
# existing revision output. It writes only to:
#   output/revision/graphon_roles_fixedK13/
# =============================================================================

options(stringsAsFactors = FALSE)

find_project_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/", mustWork = TRUE)))
  }
  if (file.exists("Mob.Rproj") || file.exists("Mob.Rmd")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  stop("Cannot locate the Mob project root. Run this script from /Mob.")
}

ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(ROOT)
SEED <- 20260811L
set.seed(SEED)

log_step <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

required_files <- c(
  "calcs/graphlists.RData",
  "graphon_distance_directed.R",
  "data/nodes.txt"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

required_packages <- c("igraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(igraph))
# Inlined graphon/geography functions are already available.
load("calcs/graphlists.RData")

if (!exists("gAGGR") || !is.list(gAGGR)) {
  stop("gAGGR was not found in calcs/graphlists.RData.")
}
if (is.null(names(gAGGR))) {
  stop("gAGGR must have year names.")
}

K_FIXED <- 13L
DESTINATION_CORE_N <- 20L
SELECTED_YEARS <- c("2015", "2018", "2019", "2020", "2021", "2024")
TRANSITIONS <- list(
  c("2015", "2018"),
  c("2018", "2019"),
  c("2019", "2020"),
  c("2020", "2021"),
  c("2021", "2024")
)

missing_years <- setdiff(SELECTED_YEARS, names(gAGGR))
if (length(missing_years)) {
  stop("Selected years missing from gAGGR: ", paste(missing_years, collapse = ", "))
}

OUT <- "output/revision/graphon_roles_fixedK13"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

nodes <- read.table(
  "data/nodes.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("nodeID", "nodeLabel", "nodeLat", "nodeLong") %in% names(nodes))) {
  stop("data/nodes.txt does not contain the expected node columns.")
}
nodes$nodeID <- as.character(nodes$nodeID)

rank01 <- function(x) {
  if (!length(x)) return(numeric())
  if (length(x) == 1L) return(0.5)
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

dense_rank_desc <- function(x) {
  ux <- sort(unique(x), decreasing = TRUE)
  match(x, ux)
}

safe_share <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  x / s
}

entropy_effective_n <- function(counts) {
  counts <- counts[counts > 0]
  if (!length(counts)) return(NA_real_)
  p <- counts / sum(counts)
  exp(-sum(p * log(p)))
}

same_block_rate <- function(labels) {
  labels <- labels[!is.na(labels)]
  n <- length(labels)
  if (n < 2L) return(NA_real_)
  M <- outer(labels, labels, FUN = "==")
  mean(M[upper.tri(M)])
}

mean_pairwise_abs <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  D <- abs(outer(x, x, FUN = "-"))
  mean(D[upper.tri(D)])
}

adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)

  tab <- table(x, y)
  comb2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))
  total_pairs <- comb2(sum(tab))
  if (total_pairs <= 0) return(NA_real_)

  expected <- sum_ai * sum_bj / total_pairs
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (abs(denom) < .Machine$double.eps) return(1)
  (sum_nij - expected) / denom
}

extract_node_metadata <- function(g) {
  ids <- igraph::V(g)$name
  if (is.null(ids)) ids <- as.character(seq_len(igraph::vcount(g)))
  ids <- as.character(ids)
  idx <- match(ids, nodes$nodeID)

  if (anyNA(idx)) {
    vlabel <- if ("nodeLabel" %in% igraph::vertex_attr_names(g)) {
      as.character(igraph::vertex_attr(g, "nodeLabel"))
    } else {
      rep(NA_character_, length(ids))
    }
    vlat <- if ("nodeLat" %in% igraph::vertex_attr_names(g)) {
      suppressWarnings(as.numeric(igraph::vertex_attr(g, "nodeLat")))
    } else {
      rep(NA_real_, length(ids))
    }
    vlong <- if ("nodeLong" %in% igraph::vertex_attr_names(g)) {
      suppressWarnings(as.numeric(igraph::vertex_attr(g, "nodeLong")))
    } else {
      rep(NA_real_, length(ids))
    }
    return(data.frame(
      nodeID = ids,
      nodeLabel = vlabel,
      nodeLat = vlat,
      nodeLong = vlong,
      stringsAsFactors = FALSE
    ))
  }

  nodes[idx, c("nodeID", "nodeLabel", "nodeLat", "nodeLong"), drop = FALSE]
}

# -----------------------------------------------------------------------------
# 1. Define a time-invariant destination-oriented reference core
# -----------------------------------------------------------------------------
log_step("Defining destination-oriented reference core from full 2006-2024 AGGR series ...")

all_years <- names(gAGGR)
reference_meta <- extract_node_metadata(gAGGR[[1]])
reference_ids <- reference_meta$nodeID

incoming_matrix <- matrix(
  NA_real_,
  nrow = length(reference_ids),
  ncol = length(all_years),
  dimnames = list(reference_ids, all_years)
)

for (year in all_years) {
  g <- gAGGR[[year]]
  A <- as_weight_matrix(g)
  meta <- extract_node_metadata(g)
  if (nrow(A) != nrow(meta)) stop("Node alignment failure while defining destination core in year ", year, ".")
  incoming_matrix[meta$nodeID, year] <- colSums(A)
}

if (anyNA(incoming_matrix)) {
  stop("Could not align all node IDs across the full AGGR time series.")
}

core_definition <- reference_meta
core_definition$mean_in_strength_2006_2024 <- rowMeans(incoming_matrix[core_definition$nodeID, , drop = FALSE])
core_definition$median_in_strength_2006_2024 <- apply(
  incoming_matrix[core_definition$nodeID, , drop = FALSE], 1, stats::median
)
core_definition$years_positive_in_strength <- rowSums(
  incoming_matrix[core_definition$nodeID, , drop = FALSE] > 0
)

ord <- order(-core_definition$mean_in_strength_2006_2024, core_definition$nodeID)
core_definition$destination_core_rank <- NA_integer_
core_definition$destination_core_rank[ord] <- seq_along(ord)
core_definition$is_destination_core <- core_definition$destination_core_rank <= DESTINATION_CORE_N
core_definition$is_regional_destination_core <-
  core_definition$is_destination_core & core_definition$nodeLabel != "Budapesti"

core_definition <- core_definition[order(core_definition$destination_core_rank), ]
rownames(core_definition) <- NULL
write.csv(
  core_definition,
  file.path(OUT, "ROLE13_destination_core_definition.csv"),
  row.names = FALSE
)

selected_core_ids <- core_definition$nodeID[core_definition$is_destination_core]
selected_regional_core_ids <- core_definition$nodeID[core_definition$is_regional_destination_core]

log_step(
  "Destination core defined: top ", DESTINATION_CORE_N,
  " regions by mean incoming weight across ", length(all_years), " years."
)

# -----------------------------------------------------------------------------
# 2. Fixed-K=13 graphon/block estimates for selected years
# -----------------------------------------------------------------------------
annual_graphons <- list()
annual_nodes <- list()
annual_blocks <- list()
annual_flows <- list()
annual_core_blocks <- list()
year_summary_rows <- list()

for (year in SELECTED_YEARS) {
  log_step("Estimating fixed K=13 block graphon for ", year, " ...")

  g <- gAGGR[[year]]
  A <- as_weight_matrix(g)
  meta <- extract_node_metadata(g)

  distinct_profiles <- n_distinct_rows(block_features(A))
  if (distinct_profiles < K_FIXED) {
    stop(
      "Year ", year, ": fixed K=", K_FIXED,
      " is infeasible because only ", distinct_profiles,
      " distinct block-feature profiles are available."
    )
  }

  go <- estimate_directed_graphon(g, K = K_FIXED, method = "block")
  if (go$K != K_FIXED) stop("Unexpected K in year ", year, ".")

  if (nrow(meta) != nrow(A) || length(go$node_labels) != nrow(A)) {
    stop("Node alignment failure in year ", year, ".")
  }

  labels <- as.integer(go$node_labels)
  sizes <- tabulate(labels, nbins = K_FIXED)
  W <- go$W
  fitted <- go$fitted

  raw_out_strength <- rowSums(A)
  raw_in_strength <- colSums(A)
  raw_out_degree <- rowSums(A > 0)
  raw_in_degree <- colSums(A > 0)
  fitted_out_mean <- rowMeans(fitted)
  fitted_in_mean <- colMeans(fitted)

  block_fitted_out_mean <- as.numeric(W %*% sizes) / sum(sizes)
  block_fitted_in_mean <- as.numeric(t(W) %*% sizes) / sum(sizes)
  block_destination_rank <- dense_rank_desc(block_fitted_in_mean)
  block_origin_rank <- dense_rank_desc(block_fitted_out_mean)
  block_destination_rank01 <- if (K_FIXED > 1L) {
    1 - (block_destination_rank - 1) / (K_FIXED - 1)
  } else {
    rep(1, K_FIXED)
  }

  core_idx <- match(meta$nodeID, core_definition$nodeID)
  if (anyNA(core_idx)) stop("Destination-core alignment failure in year ", year, ".")

  node_df <- data.frame(
    year = year,
    meta,
    block = labels,
    block_size = sizes[labels],
    block_destination_rank = block_destination_rank[labels],
    block_destination_rank01 = block_destination_rank01[labels],
    block_origin_rank = block_origin_rank[labels],
    raw_out_strength = raw_out_strength,
    raw_in_strength = raw_in_strength,
    raw_out_degree = raw_out_degree,
    raw_in_degree = raw_in_degree,
    fitted_out_mean = fitted_out_mean,
    fitted_in_mean = fitted_in_mean,
    fitted_out_rank01 = rank01(fitted_out_mean),
    fitted_in_rank01 = rank01(fitted_in_mean),
    fitted_role_balance = fitted_in_mean - fitted_out_mean,
    destination_core_rank = core_definition$destination_core_rank[core_idx],
    is_destination_core = core_definition$is_destination_core[core_idx],
    is_regional_destination_core = core_definition$is_regional_destination_core[core_idx],
    stringsAsFactors = FALSE
  )

  node_df <- node_df[order(
    node_df$block_destination_rank,
    -node_df$fitted_in_mean,
    node_df$nodeID
  ), ]
  rownames(node_df) <- NULL

  top_out_block <- vapply(seq_len(K_FIXED), function(i) which.max(W[i, ] * sizes), integer(1))
  top_in_block <- vapply(seq_len(K_FIXED), function(j) which.max(W[, j] * sizes), integer(1))
  top_out_share <- vapply(seq_len(K_FIXED), function(i) max(safe_share(W[i, ] * sizes)), numeric(1))
  top_in_share <- vapply(seq_len(K_FIXED), function(j) max(safe_share(W[, j] * sizes)), numeric(1))

  block_df <- data.frame(
    year = year,
    K = K_FIXED,
    block = seq_len(K_FIXED),
    block_size = sizes,
    block_share = sizes / sum(sizes),
    within_block_mean_weight = diag(W),
    fitted_out_mean = block_fitted_out_mean,
    fitted_in_mean = block_fitted_in_mean,
    fitted_role_balance = block_fitted_in_mean - block_fitted_out_mean,
    destination_rank = block_destination_rank,
    destination_rank01 = block_destination_rank01,
    origin_rank = block_origin_rank,
    top_out_destination_block = top_out_block,
    top_out_destination_share = top_out_share,
    top_in_origin_block = top_in_block,
    top_in_origin_share = top_in_share,
    stringsAsFactors = FALSE
  )
  block_df <- block_df[order(block_df$destination_rank, block_df$block), ]
  rownames(block_df) <- NULL

  flow_grid <- expand.grid(
    from_block = seq_len(K_FIXED),
    to_block = seq_len(K_FIXED),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  flow_grid$year <- year
  flow_grid$from_block_size <- sizes[flow_grid$from_block]
  flow_grid$to_block_size <- sizes[flow_grid$to_block]
  flow_grid$mean_weight <- W[cbind(flow_grid$from_block, flow_grid$to_block)]
  flow_grid$fitted_mass <- flow_grid$mean_weight * flow_grid$from_block_size * flow_grid$to_block_size
  total_fitted_mass <- sum(flow_grid$fitted_mass)
  flow_grid$fitted_mass_share <- if (total_fitted_mass > 0) {
    flow_grid$fitted_mass / total_fitted_mass
  } else {
    0
  }
  flow_grid <- flow_grid[order(-flow_grid$fitted_mass, flow_grid$from_block, flow_grid$to_block), ]
  rownames(flow_grid) <- NULL

  core_block_rows <- lapply(seq_len(K_FIXED), function(b) {
    in_block <- node_df$block == b
    core_nodes <- node_df[in_block & node_df$is_destination_core, , drop = FALSE]
    regional_nodes <- node_df[in_block & node_df$is_regional_destination_core, , drop = FALSE]
    data.frame(
      year = year,
      block = b,
      block_size = sizes[b],
      destination_rank = block_destination_rank[b],
      fitted_in_mean = block_fitted_in_mean[b],
      fitted_out_mean = block_fitted_out_mean[b],
      destination_core_count = nrow(core_nodes),
      destination_core_share_of_block = if (sizes[b] > 0) nrow(core_nodes) / sizes[b] else NA_real_,
      destination_core_labels = paste(core_nodes$nodeLabel, collapse = "; "),
      regional_core_count = nrow(regional_nodes),
      regional_core_share_of_block = if (sizes[b] > 0) nrow(regional_nodes) / sizes[b] else NA_real_,
      regional_core_labels = paste(regional_nodes$nodeLabel, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  core_block_df <- do.call(rbind, core_block_rows)
  core_block_df <- core_block_df[order(core_block_df$destination_rank, core_block_df$block), ]
  rownames(core_block_df) <- NULL

  core_nodes <- node_df[node_df$is_destination_core, , drop = FALSE]
  regional_nodes <- node_df[node_df$is_regional_destination_core, , drop = FALSE]
  budapest <- node_df[node_df$nodeLabel == "Budapesti", , drop = FALSE]

  core_counts <- tabulate(core_nodes$block, nbins = K_FIXED)
  regional_counts <- tabulate(regional_nodes$block, nbins = K_FIXED)

  year_summary_rows[[year]] <- data.frame(
    year = year,
    fixed_K = K_FIXED,
    n_nodes = nrow(A),
    distinct_block_feature_profiles = distinct_profiles,
    min_block_size = min(sizes),
    median_block_size = stats::median(sizes),
    max_block_size = max(sizes),
    largest_block_share = max(sizes) / sum(sizes),
    effective_number_of_all_blocks = entropy_effective_n(sizes),
    destination_core_n = nrow(core_nodes),
    destination_core_distinct_blocks = sum(core_counts > 0),
    destination_core_effective_blocks = entropy_effective_n(core_counts),
    destination_core_same_block_rate = same_block_rate(core_nodes$block),
    destination_core_pairwise_in_rank_gap = mean_pairwise_abs(core_nodes$fitted_in_rank01),
    destination_core_in_rank_sd = stats::sd(core_nodes$fitted_in_rank01),
    regional_core_n = nrow(regional_nodes),
    regional_core_distinct_blocks = sum(regional_counts > 0),
    regional_core_effective_blocks = entropy_effective_n(regional_counts),
    regional_core_same_block_rate = same_block_rate(regional_nodes$block),
    regional_core_pairwise_in_rank_gap = mean_pairwise_abs(regional_nodes$fitted_in_rank01),
    regional_core_in_rank_sd = stats::sd(regional_nodes$fitted_in_rank01),
    budapest_block_destination_rank = if (nrow(budapest) == 1L) budapest$block_destination_rank else NA_integer_,
    budapest_fitted_in_rank01 = if (nrow(budapest) == 1L) budapest$fitted_in_rank01 else NA_real_,
    raw_total_weight = sum(A),
    fitted_total_weight = sum(fitted),
    stringsAsFactors = FALSE
  )

  write.csv(node_df, file.path(OUT, paste0("ROLE13_membership_", year, ".csv")), row.names = FALSE)
  write.csv(block_df, file.path(OUT, paste0("ROLE13_block_profiles_", year, ".csv")), row.names = FALSE)
  write.csv(flow_grid, file.path(OUT, paste0("ROLE13_W_flows_", year, ".csv")), row.names = FALSE)
  write.csv(core_block_df, file.path(OUT, paste0("ROLE13_core_block_composition_", year, ".csv")), row.names = FALSE)

  annual_graphons[[year]] <- go
  annual_nodes[[year]] <- node_df
  annual_blocks[[year]] <- block_df
  annual_flows[[year]] <- flow_grid
  annual_core_blocks[[year]] <- core_block_df
}

year_summary <- do.call(rbind, year_summary_rows)
rownames(year_summary) <- NULL
write.csv(
  year_summary,
  file.path(OUT, "ROLE13_differentiation_summary_by_year.csv"),
  row.names = FALSE
)

all_nodes_selected <- do.call(rbind, annual_nodes)
rownames(all_nodes_selected) <- NULL
write.csv(
  all_nodes_selected,
  file.path(OUT, "ROLE13_memberships_selected_years_all.csv"),
  row.names = FALSE
)

all_blocks_selected <- do.call(rbind, annual_blocks)
rownames(all_blocks_selected) <- NULL
write.csv(
  all_blocks_selected,
  file.path(OUT, "ROLE13_block_profiles_selected_years_all.csv"),
  row.names = FALSE
)

core_trajectories <- all_nodes_selected[all_nodes_selected$is_destination_core, ]
core_trajectories <- core_trajectories[order(core_trajectories$destination_core_rank, core_trajectories$year), ]
rownames(core_trajectories) <- NULL
write.csv(
  core_trajectories,
  file.path(OUT, "ROLE13_destination_core_trajectories.csv"),
  row.names = FALSE
)

regional_core_trajectories <- all_nodes_selected[all_nodes_selected$is_regional_destination_core, ]
regional_core_trajectories <- regional_core_trajectories[
  order(regional_core_trajectories$destination_core_rank, regional_core_trajectories$year),
]
rownames(regional_core_trajectories) <- NULL
write.csv(
  regional_core_trajectories,
  file.path(OUT, "ROLE13_regional_core_trajectories.csv"),
  row.names = FALSE
)

all_core_blocks <- do.call(rbind, annual_core_blocks)
rownames(all_core_blocks) <- NULL
write.csv(
  all_core_blocks,
  file.path(OUT, "ROLE13_core_block_composition_all_years.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Cross-year validation without assuming block-label identity
# -----------------------------------------------------------------------------
transition_summaries <- list()
transition_node_rows <- list()
transition_overlap_rows <- list()

for (tr in TRANSITIONS) {
  from_year <- tr[1]
  to_year <- tr[2]
  log_step("Comparing fixed-K roles: ", from_year, " -> ", to_year, " ...")

  a <- annual_nodes[[from_year]]
  b <- annual_nodes[[to_year]]
  ia <- match(nodes$nodeID, a$nodeID)
  ib <- match(nodes$nodeID, b$nodeID)
  if (anyNA(ia) || anyNA(ib)) stop("Cross-year node alignment failure: ", from_year, " -> ", to_year, ".")

  a0 <- a[ia, ]
  b0 <- b[ib, ]

  ari_all <- adjusted_rand_index(a0$block, b0$block)
  core_mask <- a0$is_destination_core & b0$is_destination_core
  regional_mask <- a0$is_regional_destination_core & b0$is_regional_destination_core
  ari_core <- adjusted_rand_index(a0$block[core_mask], b0$block[core_mask])
  ari_regional <- adjusted_rand_index(a0$block[regional_mask], b0$block[regional_mask])

  node_tr <- data.frame(
    nodeID = a0$nodeID,
    nodeLabel = a0$nodeLabel,
    is_destination_core = a0$is_destination_core,
    is_regional_destination_core = a0$is_regional_destination_core,
    destination_core_rank = a0$destination_core_rank,
    from_year = from_year,
    to_year = to_year,
    block_from = a0$block,
    block_to = b0$block,
    block_destination_rank_from = a0$block_destination_rank,
    block_destination_rank_to = b0$block_destination_rank,
    fitted_out_mean_from = a0$fitted_out_mean,
    fitted_out_mean_to = b0$fitted_out_mean,
    fitted_in_mean_from = a0$fitted_in_mean,
    fitted_in_mean_to = b0$fitted_in_mean,
    fitted_out_rank01_from = a0$fitted_out_rank01,
    fitted_out_rank01_to = b0$fitted_out_rank01,
    fitted_in_rank01_from = a0$fitted_in_rank01,
    fitted_in_rank01_to = b0$fitted_in_rank01,
    raw_out_strength_from = a0$raw_out_strength,
    raw_out_strength_to = b0$raw_out_strength,
    raw_in_strength_from = a0$raw_in_strength,
    raw_in_strength_to = b0$raw_in_strength,
    stringsAsFactors = FALSE
  )

  node_tr$delta_block_destination_rank <-
    node_tr$block_destination_rank_to - node_tr$block_destination_rank_from
  node_tr$delta_fitted_out_mean <- node_tr$fitted_out_mean_to - node_tr$fitted_out_mean_from
  node_tr$delta_fitted_in_mean <- node_tr$fitted_in_mean_to - node_tr$fitted_in_mean_from
  node_tr$delta_fitted_out_rank01 <- node_tr$fitted_out_rank01_to - node_tr$fitted_out_rank01_from
  node_tr$delta_fitted_in_rank01 <- node_tr$fitted_in_rank01_to - node_tr$fitted_in_rank01_from
  node_tr$role_shift_rank_l2 <- sqrt(
    node_tr$delta_fitted_out_rank01^2 + node_tr$delta_fitted_in_rank01^2
  )

  overlap <- as.data.frame(
    table(
      from_block = factor(a0$block, levels = seq_len(K_FIXED)),
      to_block = factor(b0$block, levels = seq_len(K_FIXED))
    ),
    stringsAsFactors = FALSE
  )
  names(overlap)[3] <- "n_nodes"
  overlap$from_year <- from_year
  overlap$to_year <- to_year
  overlap$from_block <- as.integer(as.character(overlap$from_block))
  overlap$to_block <- as.integer(as.character(overlap$to_block))
  from_totals <- ave(overlap$n_nodes, overlap$from_block, FUN = sum)
  to_totals <- ave(overlap$n_nodes, overlap$to_block, FUN = sum)
  overlap$share_of_from_block <- ifelse(from_totals > 0, overlap$n_nodes / from_totals, 0)
  overlap$share_of_to_block <- ifelse(to_totals > 0, overlap$n_nodes / to_totals, 0)
  overlap <- overlap[order(-overlap$n_nodes, overlap$from_block, overlap$to_block), ]
  rownames(overlap) <- NULL

  y_from <- year_summary[year_summary$year == from_year, , drop = FALSE]
  y_to <- year_summary[year_summary$year == to_year, , drop = FALSE]

  transition_summaries[[paste0(from_year, "_", to_year)]] <- data.frame(
    from_year = from_year,
    to_year = to_year,
    fixed_K = K_FIXED,
    ARI_all_nodes = ari_all,
    ARI_destination_core = ari_core,
    ARI_regional_core = ari_regional,
    mean_abs_fitted_in_rank_shift_all = mean(abs(node_tr$delta_fitted_in_rank01), na.rm = TRUE),
    mean_abs_fitted_in_rank_shift_core = mean(abs(node_tr$delta_fitted_in_rank01[node_tr$is_destination_core]), na.rm = TRUE),
    mean_abs_fitted_in_rank_shift_regional = mean(abs(node_tr$delta_fitted_in_rank01[node_tr$is_regional_destination_core]), na.rm = TRUE),
    regional_core_effective_blocks_from = y_from$regional_core_effective_blocks,
    regional_core_effective_blocks_to = y_to$regional_core_effective_blocks,
    regional_core_same_block_rate_from = y_from$regional_core_same_block_rate,
    regional_core_same_block_rate_to = y_to$regional_core_same_block_rate,
    regional_core_pairwise_in_rank_gap_from = y_from$regional_core_pairwise_in_rank_gap,
    regional_core_pairwise_in_rank_gap_to = y_to$regional_core_pairwise_in_rank_gap,
    stringsAsFactors = FALSE
  )

  pair_name <- paste0(from_year, "_", to_year)
  write.csv(
    node_tr[order(-node_tr$role_shift_rank_l2, node_tr$destination_core_rank, node_tr$nodeID), ],
    file.path(OUT, paste0("ROLE13_node_changes_", pair_name, ".csv")),
    row.names = FALSE
  )
  write.csv(
    overlap,
    file.path(OUT, paste0("ROLE13_block_overlap_", pair_name, ".csv")),
    row.names = FALSE
  )

  transition_node_rows[[pair_name]] <- node_tr
  transition_overlap_rows[[pair_name]] <- overlap
}

transition_summary <- do.call(rbind, transition_summaries)
rownames(transition_summary) <- NULL
write.csv(
  transition_summary,
  file.path(OUT, "ROLE13_transition_summary.csv"),
  row.names = FALSE
)

all_node_changes <- do.call(rbind, transition_node_rows)
rownames(all_node_changes) <- NULL
write.csv(
  all_node_changes,
  file.path(OUT, "ROLE13_node_changes_all_transitions.csv"),
  row.names = FALSE
)

all_overlaps <- do.call(rbind, transition_overlap_rows)
rownames(all_overlaps) <- NULL
write.csv(
  all_overlaps,
  file.path(OUT, "ROLE13_block_overlap_all_transitions.csv"),
  row.names = FALSE
)

# All-pair ARI among selected years: compact label-invariant partition stability.
ari_rows <- list()
kk <- 1L
for (i in seq_len(length(SELECTED_YEARS) - 1L)) {
  for (j in seq.int(i + 1L, length(SELECTED_YEARS))) {
    y1 <- SELECTED_YEARS[i]
    y2 <- SELECTED_YEARS[j]
    a <- annual_nodes[[y1]]
    b <- annual_nodes[[y2]]
    m <- merge(
      a[, c("nodeID", "block", "is_destination_core", "is_regional_destination_core")],
      b[, c("nodeID", "block")],
      by = "nodeID",
      suffixes = c("_1", "_2"),
      sort = FALSE
    )
    ari_rows[[kk]] <- data.frame(
      year_1 = y1,
      year_2 = y2,
      ARI_all_nodes = adjusted_rand_index(m$block_1, m$block_2),
      ARI_destination_core = adjusted_rand_index(
        m$block_1[m$is_destination_core],
        m$block_2[m$is_destination_core]
      ),
      ARI_regional_core = adjusted_rand_index(
        m$block_1[m$is_regional_destination_core],
        m$block_2[m$is_regional_destination_core]
      ),
      stringsAsFactors = FALSE
    )
    kk <- kk + 1L
  }
}
ari_all_pairs <- do.call(rbind, ari_rows)
rownames(ari_all_pairs) <- NULL
write.csv(
  ari_all_pairs,
  file.path(OUT, "ROLE13_ARI_all_selected_pairs.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Output manifest and reproducibility bundle
# -----------------------------------------------------------------------------
manifest <- data.frame(
  output = c(
    "ROLE13_destination_core_definition.csv",
    "ROLE13_differentiation_summary_by_year.csv",
    "ROLE13_membership_<year>.csv",
    "ROLE13_block_profiles_<year>.csv",
    "ROLE13_W_flows_<year>.csv",
    "ROLE13_core_block_composition_<year>.csv",
    "ROLE13_destination_core_trajectories.csv",
    "ROLE13_regional_core_trajectories.csv",
    "ROLE13_transition_summary.csv",
    "ROLE13_node_changes_<from>_<to>.csv",
    "ROLE13_block_overlap_<from>_<to>.csv",
    "ROLE13_ARI_all_selected_pairs.csv"
  ),
  interpretation = c(
    "Time-invariant top-20 destination reference set defined by mean incoming application weight over 2006-2024.",
    "Primary validation table: fixed-K partition balance and differentiation metrics for all/core/regional-core nodes.",
    "Node membership and fitted/raw role measures for each selected year; block labels are snapshot-specific.",
    "Within-year block intensity profiles with label-invariant destination and origin ranks.",
    "Long-format fixed-K W matrix with fitted block-pair mass shares.",
    "Which destination-core and regional-core regions occupy each within-year structural block.",
    "Longitudinal role trajectories for the top-20 destination reference regions.",
    "Same as above excluding Budapest, for testing differentiation among regional destination centers.",
    "Label-invariant transition summary using ARI and node-aligned fitted-role rank changes.",
    "Node-level role changes across selected transitions, ranked by fitted-role displacement.",
    "Node-overlap table between snapshot-specific block partitions; does not assume label identity.",
    "Pairwise Adjusted Rand Index across every pair of selected fixed-K partitions."
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT, "ROLE13_output_manifest.csv"), row.names = FALSE)

save(
  K_FIXED,
  DESTINATION_CORE_N,
  SELECTED_YEARS,
  TRANSITIONS,
  core_definition,
  incoming_matrix,
  annual_graphons,
  annual_nodes,
  annual_blocks,
  annual_flows,
  annual_core_blocks,
  year_summary,
  transition_summary,
  ari_all_pairs,
  file = file.path(OUT, "ROLE13_validation_objects.RData")
)

capture.output(sessionInfo(), file = file.path(OUT, "ROLE13_sessionInfo.txt"))

log_step("Fixed-K=13 structural-role validation completed.")
log_step("Selected years: ", paste(SELECTED_YEARS, collapse = ", "))
log_step("Destination reference core: top ", DESTINATION_CORE_N, " by mean incoming weight, 2006-2024.")
log_step("Outputs written to ", OUT, "/")
log_step("Mob.Rmd and all existing main/revision outputs were not modified.")


})

local({
# =============================================================================
# run_graphon_multiresolution_roles.R
# Multi-resolution structural-role diagnostic for the Mob revision
# =============================================================================
#
# Scientific purpose
#   - Examine the SAME annual application network at several fixed resolutions
#     K = 2, 4, 8, 13.
#   - Test whether increasing K reveals a coarse-to-fine differentiation of
#     destination-oriented structural roles (for example, Budapest vs. the rest
#     at coarse resolution, then regional university roles at finer resolution).
#   - Do NOT interpret K itself as temporal complexity and do NOT treat graphon
#     blocks as modularity communities.
#   - Quantify how closely higher-resolution blocks refine lower-resolution
#     partitions instead of assuming hierarchical nesting.
#
# This script does NOT modify Mob.Rmd or any existing calculation objects.
# It writes only to:
#   output/revision/graphon_multiresolution/
# =============================================================================

options(stringsAsFactors = FALSE)

find_project_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/", mustWork = TRUE)))
  }
  if (file.exists("Mob.Rproj") || file.exists("Mob.Rmd")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  stop("Cannot locate the Mob project root. Run this script from /Mob.")
}

ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(ROOT)
SEED <- 20260815L
set.seed(SEED)

log_step <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

K_GRID <- c(2L, 4L, 8L, 13L)
SELECTED_YEARS <- c("2015", "2018", "2021", "2024")
MAIN_MAP_YEAR <- "2018"
OUT <- "output/revision/graphon_multiresolution"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "calcs/graphlists.RData",
  "graphon_distance_directed.R",
  "data/nodes.txt",
  "output/revision/graphon_roles_fixedK13/ROLE13_destination_core_definition.csv"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s): ", paste(missing_files, collapse = ", "))
}

required_packages <- c("igraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(igraph))
# Inlined graphon/geography functions are already available.
load("calcs/graphlists.RData")

if (!exists("gAGGR") || !is.list(gAGGR)) stop("gAGGR not found in calcs/graphlists.RData.")
if (!all(SELECTED_YEARS %in% names(gAGGR))) {
  stop("Missing selected year(s): ", paste(setdiff(SELECTED_YEARS, names(gAGGR)), collapse = ", "))
}

nodes <- read.table(
  "data/nodes.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
nodes$nodeID <- as.character(nodes$nodeID)

core_def <- read.csv(
  "output/revision/graphon_roles_fixedK13/ROLE13_destination_core_definition.csv",
  stringsAsFactors = FALSE
)
core_def$nodeID <- as.character(core_def$nodeID)

rank01 <- function(x) {
  if (!length(x)) return(numeric())
  if (length(x) == 1L) return(0.5)
  (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}

dense_rank_desc <- function(x) {
  ux <- sort(unique(x), decreasing = TRUE)
  match(x, ux)
}

entropy_effective_n <- function(labels) {
  tab <- table(labels)
  tab <- tab[tab > 0]
  if (!length(tab)) return(NA_real_)
  p <- as.numeric(tab) / sum(tab)
  exp(-sum(p * log(p)))
}

same_block_rate <- function(labels) {
  labels <- labels[!is.na(labels)]
  n <- length(labels)
  if (n < 2L) return(NA_real_)
  M <- outer(labels, labels, FUN = "==")
  mean(M[upper.tri(M)])
}

extract_node_metadata <- function(g) {
  ids <- igraph::V(g)$name
  if (is.null(ids)) ids <- as.character(seq_len(igraph::vcount(g)))
  ids <- as.character(ids)
  idx <- match(ids, nodes$nodeID)
  if (anyNA(idx)) stop("Node-ID alignment failure.")
  nodes[idx, c("nodeID", "nodeLabel", "nodeLat", "nodeLong"), drop = FALSE]
}

# Higher-resolution blocks need not be mathematically nested inside lower-K
# blocks. This score measures empirical refinement consistency:
# for each higher-K block, find its dominant parent lower-K block and compute
# the fraction of nodes belonging to that parent; the weighted overall score is
# the share of all nodes assigned consistently with such dominant parents.
refinement_diagnostic <- function(low_labels, high_labels, low_K, high_K, year) {
  tab <- table(low = low_labels, high = high_labels)
  high_sizes <- colSums(tab)
  dominant_n <- apply(tab, 2, max)
  dominant_parent <- apply(tab, 2, function(z) rownames(tab)[which.max(z)])
  high_share <- ifelse(high_sizes > 0, dominant_n / high_sizes, NA_real_)

  detail <- data.frame(
    year = year,
    low_K = low_K,
    high_K = high_K,
    high_role_rank = as.integer(colnames(tab)),
    dominant_low_role_rank = as.integer(dominant_parent),
    high_block_size = as.integer(high_sizes),
    dominant_parent_n = as.integer(dominant_n),
    dominant_parent_share = as.numeric(high_share),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    year = year,
    low_K = low_K,
    high_K = high_K,
    weighted_refinement_consistency = sum(dominant_n) / sum(tab),
    minimum_high_block_parent_share = min(high_share, na.rm = TRUE),
    median_high_block_parent_share = stats::median(high_share, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(table = tab, detail = detail, summary = summary)
}

membership_rows <- list()
block_rows <- list()
core_summary_rows <- list()
refinement_detail_rows <- list()
refinement_summary_rows <- list()
role_labels_by_year <- list()

for (year in SELECTED_YEARS) {
  log_step("Multi-resolution estimates for ", year, " ...")
  g <- gAGGR[[year]]
  A <- as_weight_matrix(g)
  meta <- extract_node_metadata(g)
  core_idx <- match(meta$nodeID, core_def$nodeID)
  if (anyNA(core_idx)) stop("Destination-core alignment failure in year ", year, ".")

  raw_out_strength <- rowSums(A)
  raw_in_strength <- colSums(A)
  year_role_labels <- list()

  for (K in K_GRID) {
    distinct_profiles <- n_distinct_rows(block_features(A))
    if (distinct_profiles < K) {
      stop("Year ", year, ": K=", K, " infeasible; only ", distinct_profiles,
           " distinct feature profiles.")
    }

    go <- estimate_directed_graphon(g, K = K, method = "block")
    labels_raw <- as.integer(go$node_labels)
    sizes_raw <- tabulate(labels_raw, nbins = K)
    W <- go$W

    # Order arbitrary k-means labels by fitted incoming intensity so that
    # role rank 1 always means the strongest destination-oriented block.
    block_fitted_in_mean_raw <- as.numeric(t(W) %*% sizes_raw) / sum(sizes_raw)
    block_fitted_out_mean_raw <- as.numeric(W %*% sizes_raw) / sum(sizes_raw)
    dest_rank_raw <- dense_rank_desc(block_fitted_in_mean_raw)
    origin_rank_raw <- dense_rank_desc(block_fitted_out_mean_raw)

    role_rank <- dest_rank_raw[labels_raw]
    origin_role_rank <- origin_rank_raw[labels_raw]
    year_role_labels[[paste0("K", K)]] <- role_rank

    fitted_in_mean <- colMeans(go$fitted)
    fitted_out_mean <- rowMeans(go$fitted)

    node_df <- data.frame(
      year = year,
      K = K,
      meta,
      raw_block = labels_raw,
      destination_role_rank = role_rank,
      origin_role_rank = origin_role_rank,
      raw_in_strength = raw_in_strength,
      raw_out_strength = raw_out_strength,
      fitted_in_mean = fitted_in_mean,
      fitted_out_mean = fitted_out_mean,
      fitted_in_rank01 = rank01(fitted_in_mean),
      destination_core_rank = core_def$destination_core_rank[core_idx],
      is_destination_core = core_def$is_destination_core[core_idx],
      is_regional_destination_core = core_def$is_regional_destination_core[core_idx],
      stringsAsFactors = FALSE
    )
    membership_rows[[length(membership_rows) + 1L]] <- node_df

    block_ids <- seq_len(K)
    block_df <- do.call(rbind, lapply(block_ids, function(b) {
      idx_b <- which(labels_raw == b)
      rr <- dest_rank_raw[b]
      in_block <- seq_len(nrow(node_df)) %in% idx_b
      core_nodes <- node_df[in_block & node_df$is_destination_core, , drop = FALSE]
      regional_nodes <- node_df[in_block & node_df$is_regional_destination_core, , drop = FALSE]
      data.frame(
        year = year,
        K = K,
        raw_block = b,
        destination_role_rank = rr,
        block_size = length(idx_b),
        fitted_in_mean = block_fitted_in_mean_raw[b],
        fitted_out_mean = block_fitted_out_mean_raw[b],
        destination_core_count = nrow(core_nodes),
        regional_core_count = nrow(regional_nodes),
        destination_core_labels = paste(core_nodes$nodeLabel, collapse = "; "),
        regional_core_labels = paste(regional_nodes$nodeLabel, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
    block_df <- block_df[order(block_df$destination_role_rank), , drop = FALSE]
    block_rows[[length(block_rows) + 1L]] <- block_df

    regional_roles <- node_df$destination_role_rank[node_df$is_regional_destination_core]
    budapest_idx <- which(node_df$nodeLabel == "Budapesti")
    core_summary_rows[[length(core_summary_rows) + 1L]] <- data.frame(
      year = year,
      K = K,
      budapest_destination_role_rank = if (length(budapest_idx)) node_df$destination_role_rank[budapest_idx[1]] else NA_integer_,
      destination_core_distinct_roles = length(unique(node_df$destination_role_rank[node_df$is_destination_core])),
      regional_core_distinct_roles = length(unique(regional_roles)),
      regional_core_effective_roles = entropy_effective_n(regional_roles),
      regional_core_same_block_rate = same_block_rate(regional_roles),
      raw_total_weight = sum(A),
      stringsAsFactors = FALSE
    )
  }

  role_labels_by_year[[year]] <- year_role_labels

  for (j in seq_len(length(K_GRID) - 1L)) {
    low_K <- K_GRID[j]
    high_K <- K_GRID[j + 1L]
    d <- refinement_diagnostic(
      year_role_labels[[paste0("K", low_K)]],
      year_role_labels[[paste0("K", high_K)]],
      low_K, high_K, year
    )
    refinement_detail_rows[[length(refinement_detail_rows) + 1L]] <- d$detail
    refinement_summary_rows[[length(refinement_summary_rows) + 1L]] <- d$summary
    write.csv(
      as.data.frame.matrix(d$table),
      file.path(OUT, paste0("MR_crosstab_", year, "_K", low_K, "_K", high_K, ".csv"))
    )
  }
}

memberships <- do.call(rbind, membership_rows)
blocks <- do.call(rbind, block_rows)
core_summary <- do.call(rbind, core_summary_rows)
refinement_detail <- do.call(rbind, refinement_detail_rows)
refinement_summary <- do.call(rbind, refinement_summary_rows)

write.csv(memberships, file.path(OUT, "MR_memberships_all.csv"), row.names = FALSE)
write.csv(blocks, file.path(OUT, "MR_block_profiles_all.csv"), row.names = FALSE)
write.csv(core_summary, file.path(OUT, "MR_core_summary.csv"), row.names = FALSE)
write.csv(refinement_detail, file.path(OUT, "MR_refinement_detail.csv"), row.names = FALSE)
write.csv(refinement_summary, file.path(OUT, "MR_refinement_summary.csv"), row.names = FALSE)

# Optional publication-oriented map output. It reuses the project's existing
# geography and geoplot helper when available. Failure here does NOT invalidate
# the numerical diagnostic above.
if (file.exists("hungary.RData") && file.exists("geoplot.R")) {
  try({
    load("hungary.RData")
# Inlined graphon/geography functions are already available.
    if (exists("shp")) {
      fixed_palette <- grDevices::hcl.colors(13, palette = "Dark 3")

      draw_multires_panel <- function(year, K) {
        g <- gAGGR[[year]]
        df <- memberships[memberships$year == year & memberships$K == K, , drop = FALSE]
        ids <- as.character(igraph::V(g)$name)
        o <- match(ids, df$nodeID)
        if (anyNA(o)) stop("Map node alignment failure for year ", year, ", K=", K)
        df <- df[o, , drop = FALSE]

        igraph::V(g)$lat <- df$nodeLat
        igraph::V(g)$lon <- df$nodeLong

        # Rank 1 is the strongest destination role. Map ranks over the full
        # 13-color scale so colors carry comparable coarse-to-fine meaning.
        color_idx <- if (K == 1L) 1L else 1L + round((df$destination_role_rank - 1) * 12 / (K - 1))
        node_colors <- fixed_palette[pmax(1L, pmin(13L, color_idx))]

        geoplot(
          g,
          node_col = node_colors,
          mode = "custom_shape",
          shape = shp,
          shape_bg = "white",
          shape_border_col = "gray70",
          node_size = df$raw_in_strength,
          node_size_scale = c(0.6, 2.4),
          edge_alpha = 0.05,
          edge_width = 0.35,
          edge_curved = TRUE,
          show_labels = TRUE,
          label_size_mult = 0.45,
          top_n_labels = 13,
          arrow_size = 0
        )
        title(main = paste0(year, ": K=", K), line = -1.1, cex.main = 0.95)
      }

      grDevices::pdf(
        file.path(OUT, paste0("MR_multiresolution_map_", MAIN_MAP_YEAR, ".pdf")),
        width = 10.5, height = 8.2, onefile = TRUE
      )
      oldpar <- par(no.readonly = TRUE)
      par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.2, 0.5))
      for (K in K_GRID) draw_multires_panel(MAIN_MAP_YEAR, K)
      par(oldpar)
      grDevices::dev.off()

      # Supplementary: the same four resolutions for all prespecified years.
      grDevices::pdf(
        file.path(OUT, "MR_multiresolution_maps_selected_years.pdf"),
        width = 10.5, height = 8.2, onefile = TRUE
      )
      for (year in SELECTED_YEARS) {
        oldpar <- par(no.readonly = TRUE)
        par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.2, 0.5))
        for (K in K_GRID) draw_multires_panel(year, K)
        par(oldpar)
      }
      grDevices::dev.off()
    }
  }, silent = TRUE)
}

save(
  memberships, blocks, core_summary, refinement_detail, refinement_summary,
  K_GRID, SELECTED_YEARS, MAIN_MAP_YEAR,
  file = file.path(OUT, "MR_multiresolution_objects.RData")
)

writeLines(capture.output(sessionInfo()), file.path(OUT, "MR_sessionInfo.txt"))

cat("\n=== MULTI-RESOLUTION RUN COMPLETED ===\n")
cat("Outputs:", OUT, "\n")


})

local({
# =============================================================================
# run_network_role_flow_figure_v2.R
# Network geography of major interregional application flows embedded in the
# fixed-K=13 graphon structural-role representation.
# =============================================================================

options(stringsAsFactors = FALSE)

find_project_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep('^--file=', a, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub('^--file=', '', f[1]), winslash = '/', mustWork = TRUE)))
  }
  if (file.exists('Mob.Rproj') || file.exists('Mob.Rmd') || file.exists('Mob_revised.Rmd') || file.exists('Mob_revised_v2.Rmd')) {
    return(normalizePath('.', winslash = '/', mustWork = TRUE))
  }
  stop('Cannot locate the Mob project root. Run this script from /Mob.')
}

ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(ROOT)

log_step <- function(...) {
  cat(sprintf('[%s] ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')), ..., '\n', sep = '')
}

YEAR_MAIN <- '2018'
OUTDIR <- 'output/revision/network_role_flows'
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  'calcs/graphlists.RData',
  file.path('output/revision/graphon_roles_fixedK13', paste0('ROLE13_membership_', YEAR_MAIN, '.csv'))
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop('Missing required file(s): ', paste(missing_files, collapse = ', '))

required_packages <- c('igraph', 'ggplot2')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop('Missing required R package(s): ', paste(missing_packages, collapse = ', '))

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggplot2))
use_maps <- requireNamespace('maps', quietly = TRUE)
use_ggrepel <- requireNamespace('ggrepel', quietly = TRUE)

normalize_label <- function(x) {
  y <- tolower(as.character(x))
  y <- iconv(y, from = '', to = 'ASCII//TRANSLIT')
  gsub('[^a-z0-9]+', '', y)
}

log_step('Loading the 2018 application network and fixed-K=13 role membership ...')
load('calcs/graphlists.RData')
if (!exists('gAGGR') || !is.list(gAGGR) || is.null(gAGGR[[YEAR_MAIN]])) {
  stop('gAGGR or the selected year was not found in calcs/graphlists.RData.')
}

roles <- read.csv(
  file.path('output/revision/graphon_roles_fixedK13', paste0('ROLE13_membership_', YEAR_MAIN, '.csv')),
  stringsAsFactors = FALSE
)
required_role_cols <- c(
  'nodeID','nodeLabel','nodeLat','nodeLong','block','block_destination_rank',
  'fitted_in_mean','fitted_out_mean','is_destination_core','is_regional_destination_core'
)
missing_role_cols <- setdiff(required_role_cols, names(roles))
if (length(missing_role_cols)) stop('Role file missing column(s): ', paste(missing_role_cols, collapse = ', '))
roles$nodeID <- as.character(roles$nodeID)
roles$nodeLabel_norm <- normalize_label(roles$nodeLabel)
roles$role_group <- paste0('Block ', roles$block)
roles$destination_tier <- paste0('Tier ', roles$block_destination_rank)

# Directed edge: applicant-residence micro-region -> chosen-campus micro-region.
g <- gAGGR[[YEAR_MAIN]]
edges <- igraph::as_data_frame(g, what = 'edges')
if (!'weight' %in% names(edges)) stop('The selected graph has no edge weight attribute.')
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)
edges <- edges[edges$from != edges$to, , drop = FALSE]

from_meta <- roles[, c('nodeID','nodeLabel','nodeLabel_norm','nodeLat','nodeLong','block','block_destination_rank','role_group','destination_tier','is_destination_core','is_regional_destination_core')]
names(from_meta) <- c('from','from_label','from_label_norm','from_lat','from_long','from_block','from_tier_rank','from_role_group','from_destination_tier','from_destination_core','from_regional_core')
to_meta <- roles[, c('nodeID','nodeLabel','nodeLabel_norm','nodeLat','nodeLong','block','block_destination_rank','role_group','destination_tier','is_destination_core','is_regional_destination_core')]
names(to_meta) <- c('to','to_label','to_label_norm','to_lat','to_long','to_block','to_tier_rank','to_role_group','to_destination_tier','to_destination_core','to_regional_core')

edges <- merge(edges, from_meta, by = 'from', all.x = TRUE, sort = FALSE)
edges <- merge(edges, to_meta, by = 'to', all.x = TRUE, sort = FALSE)
if (anyNA(edges$from_block) || anyNA(edges$to_block)) stop('Could not align all edge endpoints to fixed-K=13 roles.')

edges$same_role_group <- edges$from_block == edges$to_block
edges$same_destination_tier <- edges$from_tier_rank == edges$to_tier_rank
edges$same_tier_different_group <- edges$same_destination_tier & !edges$same_role_group
edges$cross_tier <- !edges$same_destination_tier
edges$either_destination_core <- edges$from_destination_core | edges$to_destination_core
edges$both_regional_core <- edges$from_regional_core & edges$to_regional_core
edges$pair_id <- paste(pmin(edges$from_label_norm, edges$to_label_norm), pmax(edges$from_label_norm, edges$to_label_norm), sep = '__')
edges <- edges[order(-edges$weight), ]

write.csv(edges, file.path(OUTDIR, paste0('NRF_all_interregional_edges_', YEAR_MAIN, '.csv')), row.names = FALSE)
write.csv(roles, file.path(OUTDIR, paste0('NRF_node_roles_', YEAR_MAIN, '.csv')), row.names = FALSE)

# Top links by structural relationship.
top_n <- function(d, n = 30L) head(d[order(-d$weight), , drop = FALSE], n)
within_group_top <- top_n(edges[edges$same_role_group, , drop = FALSE], 30)
same_tier_diff_top <- top_n(edges[edges$same_tier_different_group, , drop = FALSE], 30)
cross_tier_top <- top_n(edges[edges$cross_tier, , drop = FALSE], 30)
write.csv(within_group_top, file.path(OUTDIR, paste0('NRF_top_within_role_group_edges_', YEAR_MAIN, '.csv')), row.names = FALSE)
write.csv(same_tier_diff_top, file.path(OUTDIR, paste0('NRF_top_same_tier_different_group_edges_', YEAR_MAIN, '.csv')), row.names = FALSE)
write.csv(cross_tier_top, file.path(OUTDIR, paste0('NRF_top_cross_tier_edges_', YEAR_MAIN, '.csv')), row.names = FALSE)

# Selected university micro-regions for interpretive examples.
selected_city_norm <- c('gyori','veszpremi','pecsi','szegedi','debreceni','miskolci','godolloi','budapesti')
selected_flows <- edges[
  edges$from_label_norm %in% selected_city_norm & edges$to_label_norm %in% selected_city_norm,
  , drop = FALSE
]
selected_flows <- selected_flows[order(-selected_flows$weight), ]
write.csv(
  selected_flows,
  file.path(OUTDIR, paste0('NRF_selected_city_pair_flows_', YEAR_MAIN, '_v3.csv')),
  row.names = FALSE
)

# Pair-level summary, including the motivating Gyor-Veszprem example.
selected_pairs <- list(
  c('gyori','veszpremi'),
  c('pecsi','szegedi'),
  c('pecsi','debreceni'),
  c('szegedi','debreceni'),
  c('gyori','godolloi'),
  c('veszpremi','godolloi')
)
pair_rows <- lapply(selected_pairs, function(z) {
  pid <- paste(sort(z), collapse = '__')
  d <- edges[edges$pair_id == pid, , drop = FALSE]
  node_a <- roles[roles$nodeLabel_norm == z[1], , drop = FALSE]
  node_b <- roles[roles$nodeLabel_norm == z[2], , drop = FALSE]
  a_label <- if (nrow(node_a)) node_a$nodeLabel[1] else z[1]
  b_label <- if (nrow(node_b)) node_b$nodeLabel[1] else z[2]
  ab <- d[d$from_label_norm == z[1] & d$to_label_norm == z[2], , drop = FALSE]
  ba <- d[d$from_label_norm == z[2] & d$to_label_norm == z[1], , drop = FALSE]
  data.frame(
    year = YEAR_MAIN,
    pair = paste(a_label, b_label, sep = ' <-> '),
    A_to_B = if (nrow(ab)) sum(ab$weight) else 0,
    B_to_A = if (nrow(ba)) sum(ba$weight) else 0,
    bidirectional_total = sum(d$weight),
    same_role_group = if (nrow(node_a) && nrow(node_b)) node_a$block[1] == node_b$block[1] else NA,
    same_destination_tier = if (nrow(node_a) && nrow(node_b)) node_a$block_destination_rank[1] == node_b$block_destination_rank[1] else NA,
    A_role_group = if (nrow(node_a)) node_a$role_group[1] else NA,
    B_role_group = if (nrow(node_b)) node_b$role_group[1] else NA,
    A_destination_tier = if (nrow(node_a)) node_a$destination_tier[1] else NA,
    B_destination_tier = if (nrow(node_b)) node_b$destination_tier[1] else NA,
    stringsAsFactors = FALSE
  )
})
pair_summary <- do.call(rbind, pair_rows)
pair_summary <- pair_summary[order(-pair_summary$bidirectional_total), ]
write.csv(
  pair_summary,
  file.path(OUTDIR, paste0('NRF_selected_pair_summary_', YEAR_MAIN, '_v3.csv')),
  row.names = FALSE
)

# Global descriptive shares. These refer to interregional application weight only.
total_w <- sum(edges$weight)
regional_core <- roles[roles$is_regional_destination_core, , drop = FALSE]
regional_effective_roles <- NA_real_
if (nrow(regional_core)) {
  tab <- table(regional_core$block)
  p <- as.numeric(tab / sum(tab))
  regional_effective_roles <- exp(-sum(p * log(p)))
}
summary_tab <- data.frame(
  year = YEAR_MAIN,
  interregional_edges = nrow(edges),
  interregional_total_weight = total_w,
  weighted_share_same_role_group = sum(edges$weight[edges$same_role_group]) / total_w,
  weighted_share_same_tier_different_group = sum(edges$weight[edges$same_tier_different_group]) / total_w,
  weighted_share_any_same_destination_tier = sum(edges$weight[edges$same_destination_tier]) / total_w,
  weighted_share_cross_tier = sum(edges$weight[edges$cross_tier]) / total_w,
  regional_core_count = nrow(regional_core),
  regional_core_distinct_role_groups = length(unique(regional_core$block)),
  regional_core_distinct_destination_tiers = length(unique(regional_core$block_destination_rank)),
  regional_core_effective_role_groups = regional_effective_roles,
  stringsAsFactors = FALSE
)
write.csv(summary_tab, file.path(OUTDIR, paste0('NRF_summary_', YEAR_MAIN, '.csv')), row.names = FALSE)

# Plot subset: strongest general links + additional within-group/same-tier links + selected pairs.
plot_general <- top_n(edges, 120)
plot_within <- top_n(edges[edges$same_role_group, , drop = FALSE], 45)
plot_edges <- unique(rbind(plot_general, plot_within, selected_flows))
selected_pair_ids <- unique(selected_flows$pair_id)
plot_edges$edge_class <- ifelse(
  plot_edges$pair_id %in% selected_pair_ids,
  'Selected university-city pair',
  ifelse(plot_edges$same_role_group,
         'Within same structural-role group',
         'Across destination tiers')
)
plot_edges$edge_class <- factor(
  plot_edges$edge_class,
  levels = c('Across destination tiers','Within same structural-role group','Selected university-city pair')
)

hun_map <- NULL
if (use_maps) {
  hun_map <- tryCatch(ggplot2::map_data('world', region = 'Hungary'), error = function(e) NULL)
}

label_keep <- c('budapesti','debreceni','szegedi','pecsi','gyori','veszpremi','miskolci','godolloi')
label_nodes <- roles[roles$nodeLabel_norm %in% label_keep, , drop = FALSE]
roles$destination_tier_factor <- factor(
  roles$destination_tier,
  levels = paste0('Tier ', sort(unique(roles$block_destination_rank)))
)

log_step('Rendering the directed network geography figure ...')
p <- ggplot()
if (!is.null(hun_map) && nrow(hun_map)) {
  p <- p + geom_polygon(
    data = hun_map,
    aes(x = long, y = lat, group = group),
    fill = 'gray97', color = 'gray80', linewidth = 0.25
  )
}
p <- p +
  geom_segment(
    data = plot_edges,
    aes(x = from_long, y = from_lat, xend = to_long, yend = to_lat,
        linewidth = weight, alpha = weight, color = edge_class),
    arrow = grid::arrow(length = grid::unit(0.055, 'inches'), type = 'closed'),
    lineend = 'round'
  ) +
  geom_point(
    data = roles,
    aes(x = nodeLong, y = nodeLat, fill = destination_tier_factor),
    shape = 21, size = 2.8, stroke = 0.25, color = 'black'
  ) +
  scale_linewidth_continuous(range = c(0.15, 1.7), guide = 'none') +
  scale_alpha_continuous(range = c(0.15, 0.80), guide = 'none') +
  scale_color_manual(values = c(
    'Across destination tiers' = 'gray60',
    'Within same structural-role group' = '#2C7FB8',
    'Selected university-city pair' = '#D95F0E'
  )) +
  scale_fill_viridis_d(name = 'K = 13 destination-role tier', option = 'plasma', direction = -1) +
  coord_equal() +
  labs(
    title = 'Network geography of major interregional application flows (2018)',
    subtitle = 'Micro-regions are colored by fixed-K = 13 destination-role tier; only the strongest interregional application flows are displayed',
    x = NULL, y = NULL, color = 'Edge type',
    caption = 'Graphon blocks represent structural-role groups, not modularity-based communities. Edge direction: applicant residence -> chosen campus.'
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = 'right', plot.title = element_text(face = 'bold'), plot.caption = element_text(size = 8)
  )

if (nrow(label_nodes)) {
  if (use_ggrepel) {
    p <- p + ggrepel::geom_text_repel(
      data = label_nodes,
      aes(x = nodeLong, y = nodeLat, label = nodeLabel),
      size = 3, min.segment.length = 0, box.padding = 0.2, point.padding = 0.15, seed = 123
    )
  } else {
    p <- p + geom_text(
      data = label_nodes,
      aes(x = nodeLong, y = nodeLat, label = nodeLabel),
      size = 3, nudge_y = 0.10, check_overlap = TRUE
    )
  }
}

ggsave(
  file.path(OUTDIR, paste0('NRF_network_map_', YEAR_MAIN, '_v3.pdf')),
  p, width = 11, height = 7.2, units = 'in'
)
# PNG export intentionally omitted; publication outputs are vector PDF.

log_step('Completed. Outputs written to ', OUTDIR, '/')


})

# =============================================================================
# 9. Vector-PDF manuscript figures (same data and logic as Mob_final.Rmd)
# =============================================================================

save_pdf <- function(plot, filename, width, height) {
  ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                  units = "in", device = grDevices::pdf)
}

data_melted <- reshape2::melt(fig01, id.vars = "Professions",
                              variable.name = "Year", value.name = "Applicants")
p_descr <- ggplot2::ggplot(data_melted,
  ggplot2::aes(fill = Professions, y = Applicants, x = Year)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") +
  ggplot2::scale_fill_brewer(palette = "Paired") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom", legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 5),
    legend.key.height = grid::unit(0.4, "lines"),
    legend.key.width = grid::unit(0.9, "lines")
  )
save_pdf(p_descr, "descr-1.pdf", 8, 6)

zlong <- tidyr::pivot_longer(z, cols = dplyr::all_of(MAIN_NETWORK_INDICATORS),
                             names_to = "Indicator", values_to = "z")
p_ts <- ggplot2::ggplot(zlong, ggplot2::aes(x = Year, y = z, group = Indicator)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.0) +
  ggplot2::facet_wrap(~Indicator, ncol = 2) +
  ggplot2::scale_x_continuous(breaks = seq(2006, 2024, 3)) +
  ggplot2::labs(x = NULL, y = "Standardized value") +
  ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "none")
save_pdf(p_ts, "tsall-1.pdf", 8, 7)

FINAL_EST <- utils::read.csv("output/revision/FINAL_fixedK13_estimator_trajectories.csv")
FINAL_EST_AGREE <- utils::read.csv("output/revision/FINAL_fixedK13_estimator_agreement.csv")
FINAL_KSENS <- utils::read.csv("output/revision/FINAL_K_sensitivity_trajectories.csv")
FINAL_KAGREE <- utils::read.csv("output/revision/FINAL_K_sensitivity_vs_K13.csv")
ROLE13_SUMMARY <- utils::read.csv(
  "output/revision/graphon_roles_fixedK13/ROLE13_differentiation_summary_by_year.csv")
ROLE13_TRANSITIONS <- utils::read.csv(
  "output/revision/graphon_roles_fixedK13/ROLE13_transition_summary.csv")
MR_MEMBERSHIP <- utils::read.csv(
  "output/revision/graphon_multiresolution/MR_memberships_all.csv")

ks <- FINAL_KSENS |>
  dplyr::filter(metric %in% c("frobenius","cut","wasserstein"),
                fixed_K %in% c(3,4,5,6,8,10,13)) |>
  dplyr::mutate(Year = as.integer(sub(".*-", "", year_pair)))
band <- ks |>
  dplyr::group_by(metric, Year) |>
  dplyr::summarise(low = min(distance, na.rm = TRUE), high = max(distance, na.rm = TRUE),
                   median = stats::median(distance, na.rm = TRUE), .groups = "drop")
ref <- dplyr::filter(ks, fixed_K == 13)
p_grall <- ggplot2::ggplot(band, ggplot2::aes(x = Year)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high), alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = median), linetype = 2, linewidth = 0.45) +
  ggplot2::geom_line(data = ref, ggplot2::aes(y = distance), linewidth = 0.8) +
  ggplot2::geom_point(data = ref, ggplot2::aes(y = distance), size = 1.2) +
  ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y",
    labeller = ggplot2::as_labeller(c(
      frobenius = "Frobenius (primary)", cut = "Cut (primary)",
      wasserstein = "Wasserstein (secondary)"))) +
  ggplot2::scale_x_continuous(breaks = 2007:2024) +
  ggplot2::labs(x = "Later year in consecutive pair", y = "Distance") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_pdf(p_grall, "grall-1.pdf", 8, 7)

est <- FINAL_EST |>
  dplyr::filter(metric %in% c("frobenius","cut")) |>
  dplyr::mutate(
    Year = as.integer(sub(".*-", "", year_pair)),
    Estimator = factor(estimator, levels = c("block","smooth","spectral"))
  )
p_heat <- ggplot2::ggplot(est,
  ggplot2::aes(x = Year, y = distance, linetype = Estimator, group = Estimator)) +
  ggplot2::geom_line(linewidth = 0.65) + ggplot2::geom_point(size = 1.1) +
  ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y",
                      labeller = ggplot2::as_labeller(c(frobenius = "Frobenius", cut = "Cut"))) +
  ggplot2::scale_x_continuous(breaks = 2007:2024) +
  ggplot2::labs(x = "Later year in consecutive pair", y = "Distance") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_pdf(p_heat, "grheatmap-1.pdf", 8, 6)

mr18 <- MR_MEMBERSHIP |>
  dplyr::filter(year == 2018, K %in% c(2,4,8,13)) |>
  dplyr::mutate(
    K = factor(K, levels = c(2,4,8,13), labels = paste0("K=", c(2,4,8,13))),
    Role = factor(destination_role_rank)
  )
mr18_lab <- dplyr::filter(mr18, destination_core_rank <= 10)
p_mr <- ggplot2::ggplot(mr18, ggplot2::aes(x = nodeLong, y = nodeLat)) +
  ggplot2::geom_point(ggplot2::aes(color = Role), size = 1.35, alpha = 0.82) +
  ggplot2::geom_point(data = dplyr::filter(mr18, is_destination_core),
                      shape = 21, fill = NA, color = "black", size = 2.2, stroke = 0.35) +
  ggplot2::geom_text(data = mr18_lab, ggplot2::aes(label = nodeLabel),
                     size = 2.0, check_overlap = TRUE, nudge_y = 0.07) +
  ggplot2::facet_wrap(~K, ncol = 2) + ggplot2::coord_quickmap() +
  ggplot2::labs(x = NULL, y = NULL, color = "Destination role rank") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
save_pdf(p_mr, "grmultiresolution-1.pdf", 8, 7)
save_pdf(p_mr,
  "output/revision/graphon_multiresolution/MR_multiresolution_map_2018.pdf", 10.5, 8.2)

rd <- ROLE13_SUMMARY |>
  dplyr::select(year, regional_core_effective_blocks, regional_core_same_block_rate,
                regional_core_pairwise_in_rank_gap) |>
  tidyr::pivot_longer(-year, names_to = "measure", values_to = "value")
rd$measure <- factor(
  rd$measure,
  levels = c("regional_core_effective_blocks","regional_core_same_block_rate",
             "regional_core_pairwise_in_rank_gap"),
  labels = c("Effective number of occupied blocks","Same-block pair rate",
             "Mean pairwise destination-rank gap")
)
p_role <- ggplot2::ggplot(rd, ggplot2::aes(x = year, y = value)) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~measure, ncol = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = ROLE13_SUMMARY$year) +
  ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_bw(base_size = 9)
save_pdf(p_role, "gr15-1.pdf", 8, 6)

tr <- ROLE13_TRANSITIONS |>
  dplyr::mutate(Transition = paste0(from_year, "-", to_year)) |>
  dplyr::select(Transition, ARI_regional_core, mean_abs_fitted_in_rank_shift_regional) |>
  tidyr::pivot_longer(-Transition, names_to = "measure", values_to = "value")
tr$measure <- factor(
  tr$measure,
  levels = c("ARI_regional_core","mean_abs_fitted_in_rank_shift_regional"),
  labels = c("Adjusted Rand index: regional destination core",
             "Mean absolute fitted destination-rank shift")
)
p_turn <- ggplot2::ggplot(tr, ggplot2::aes(x = Transition, y = value, group = measure)) +
  ggplot2::geom_col() + ggplot2::facet_wrap(~measure, ncol = 1, scales = "free_y") +
  ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_pdf(p_turn, "grgeos-1.pdf", 8, 5.5)
log_step("Saved vector-PDF manuscript figures.")

# =============================================================================
# 10. Reproducibility metadata, scientific checks and output assertions
# =============================================================================

validation <- data.frame(
  check = c("years_2006_2024", "nodes_175", "year_2019_total_weight",
            "main_multiresolution_year_2018"),
  observed = c(
    paste(names(gAGGR), collapse = ","),
    as.character(igraph::vcount(gAGGR[[1L]])),
    as.character(sum(as.numeric(igraph::E(gAGGR[["2019"]])$weight), na.rm = TRUE)),
    "2018"
  ),
  expected = c(paste(2006:2024, collapse = ","), "175", "384035", "2018"),
  stringsAsFactors = FALSE
)
validation$passed <- validation$observed == validation$expected
utils::write.csv(validation, "calcs/run_all_standalone_validation.csv", row.names = FALSE)
if (!all(validation$passed)) {
  warning("One or more scientific invariant checks failed; inspect calcs/run_all_standalone_validation.csv.")
}

required_outputs <- c(
  "calcs/graphlists.RData", "calcs/netprops.RData", "calcs/nodeprops.RData",
  "calcs/graphon_dists.RData", "calcs/revision/FINAL_graphon_revision_objects.RData",
  "output/FIG01.xlsx", "output/FIG02_MAIN_NETWORK_INDICATORS.xlsx",
  "output/FIG03_FOCUSED_NETWORK_TRAJECTORIES.xlsx",
  "output/FIG03_MAIN_NODE_INDICATORS.xlsx", "output/FIGA1_NETIND_ALL.xlsx",
  "output/FIGA2_PREFERENCES.xlsx", "output/FIGA3_PROFESSIONS.xlsx",
  "output/revision/FINAL_graphon_specification.csv",
  "output/revision/FINAL_fixedK13_estimator_trajectories.csv",
  "output/revision/FINAL_fixedK13_estimator_agreement.csv",
  "output/revision/FINAL_K_sensitivity_trajectories.csv",
  "output/revision/FINAL_K_sensitivity_vs_K13.csv",
  "output/revision/FINAL_primary_transition_ranking.csv",
  "output/revision/FINAL_singular_spectrum_benchmark.csv",
  "output/revision/graphon_roles_fixedK13/ROLE13_differentiation_summary_by_year.csv",
  "output/revision/graphon_roles_fixedK13/ROLE13_transition_summary.csv",
  "output/revision/graphon_multiresolution/MR_core_summary.csv",
  "output/revision/graphon_multiresolution/MR_refinement_summary.csv",
  "output/revision/graphon_multiresolution/MR_block_profiles_all.csv",
  "output/revision/graphon_multiresolution/MR_memberships_all.csv",
  "output/revision/network_role_flows/NRF_all_interregional_edges_2018.csv",
  "output/revision/network_role_flows/NRF_node_roles_2018.csv",
  "output/revision/network_role_flows/NRF_top_within_role_group_edges_2018.csv",
  "output/revision/network_role_flows/NRF_top_same_tier_different_group_edges_2018.csv",
  "output/revision/network_role_flows/NRF_top_cross_tier_edges_2018.csv",
  "output/revision/network_role_flows/NRF_summary_2018.csv",
  "output/revision/network_role_flows/NRF_selected_city_pair_flows_2018_v3.csv",
  "output/revision/network_role_flows/NRF_selected_pair_summary_2018_v3.csv",
  "output/revision/network_role_flows/NRF_network_map_2018_v3.pdf",
  "descr-1.pdf", "tsall-1.pdf", "grall-1.pdf", "grheatmap-1.pdf",
  "grmultiresolution-1.pdf", "gr15-1.pdf", "grgeos-1.pdf",
  "hungary.RData", "tokable.R", "geoplot.R", "graphon_distance_directed.R"
)
assert_files(required_outputs)

input_files <- c("run_all_standalone.R", list.files("data", full.names = TRUE))
input_files <- input_files[file.exists(input_files)]
input_manifest <- data.frame(
  file = input_files,
  bytes = unname(file.info(input_files)$size),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, "calcs/run_all_standalone_input_manifest.csv", row.names = FALSE)

generated <- unique(c(
  required_outputs,
  list.files("calcs", recursive = TRUE, full.names = TRUE),
  list.files("output", recursive = TRUE, full.names = TRUE)
))
generated <- generated[file.exists(generated) & !dir.exists(generated)]
output_manifest <- data.frame(
  file = generated,
  bytes = unname(file.info(generated)$size),
  modified = as.character(file.info(generated)$mtime),
  md5 = unname(tools::md5sum(generated)),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, "calcs/run_all_standalone_output_manifest.csv", row.names = FALSE)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(required_packages, function(p) as.character(utils::packageVersion(p)), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(package_versions, "calcs/run_all_standalone_package_versions.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "calcs/run_all_standalone_sessionInfo.txt")
writeLines(c(
  paste0("global_seed=", SEED),
  "graphon_internal_kmeans_seed=42",
  "cut_distance_seed=42",
  paste0("network_lacunarity_repetitions=", NLAC_REPS),
  "multi_resolution_component_seed=20260815"
), "calcs/run_all_standalone_seed_specification.txt")

log_step("Standalone pipeline completed successfully.")
log_step("Mob_final.Rmd can now be rendered with bibliography-final.bib present.")

# Empirical benchmark: transition evidence, robustness, and role interpretation.

mx_log("Stage 03: running the empirical benchmark.")
mx_require_packages(c("igraph"))

prepared_path <- file.path(
  MX_PATHS$derived_empirical, "EMPIRICAL_node_aligned_networks.rds"
)
mx_assert_files(prepared_path)
prepared <- readRDS(prepared_path)
networks <- prepared$networks
nodes <- prepared$nodes
years <- names(networks)
adjacency <- lapply(networks, as_weight_matrix)

annual_indicators <- do.call(rbind, lapply(years, function(year) {
  values <- mx_network_indicators(adjacency[[year]])
  data.frame(
    year = as.integer(year),
    indicator = names(values),
    value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}))
mx_write_csv(
  annual_indicators,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_indicators.csv")
)

fits_path <- file.path(MX_PATHS$object_output, "EMPIRICAL_graphon_fits.rds")
if (MX_CONFIG$reuse_intermediate && file.exists(fits_path)) {
  mx_log("Reusing empirical graphon fits from ", fits_path, ".")
  fits <- readRDS(fits_path)
  expected_signature <- list(
    years = years,
    reference_K = MX_CONFIG$reference_K,
    estimators = MX_CONFIG$estimators,
    resolution_grid = MX_CONFIG$resolution_grid
  )
  if (!identical(fits$signature, expected_signature)) {
    stop(
      "Existing empirical graphon fits do not match the current configuration. ",
      "Set METHODSX_REUSE_INTERMEDIATE=0 and rerun.", call. = FALSE
    )
  }
} else {
  fixed_estimator_fits <- setNames(
    vector("list", length(MX_CONFIG$estimators)), MX_CONFIG$estimators
  )
  for (method in MX_CONFIG$estimators) {
    mx_log("Fitting empirical estimator: ", method, ", K=", MX_CONFIG$reference_K, ".")
    fixed_estimator_fits[[method]] <- setNames(lapply(years, function(year) {
      set.seed(MX_CONFIG$seed + as.integer(year))
      mx_fit_graphon(adjacency[[year]], K = MX_CONFIG$reference_K, method = method)
    }), years)
  }

  resolution_fits <- setNames(
    vector("list", length(MX_CONFIG$resolution_grid)),
    paste0("K", MX_CONFIG$resolution_grid)
  )
  for (K in MX_CONFIG$resolution_grid) {
    key <- paste0("K", K)
    if (K == MX_CONFIG$reference_K) {
      resolution_fits[[key]] <- fixed_estimator_fits$block
    } else {
      mx_log("Fitting empirical block estimator at K=", K, ".")
      resolution_fits[[key]] <- setNames(lapply(years, function(year) {
        set.seed(MX_CONFIG$seed + 1000L * K + as.integer(year))
        mx_fit_graphon(adjacency[[year]], K = K, method = "block")
      }), years)
    }
  }

  fits <- list(
    signature = list(
      years = years,
      reference_K = MX_CONFIG$reference_K,
      estimators = MX_CONFIG$estimators,
      resolution_grid = MX_CONFIG$resolution_grid
    ),
    fixed_estimator = fixed_estimator_fits,
    resolution = resolution_fits
  )
  saveRDS(fits, fits_path, compress = "gzip")
}

metric_family <- function(channel) {
  if (grepl("^indicator_", channel)) return("conventional indicator")
  if (grepl("^graphon_", channel)) return("fitted-kernel distance")
  if (grepl("^singular_", channel)) return("adjacency benchmark")
  if (channel == "estimated_role_instability") return("structural roles")
  "other"
}

metric_scale <- function(channel) {
  if (grepl("_shape$", channel)) return("unit-mean shape")
  if (grepl("_raw$", channel)) return("raw weighted scale")
  if (channel == "estimated_role_instability") return("label-invariant partition")
  "native indicator scale"
}

transition_rows <- list()
row_id <- 1L
append_transition <- function(from_year, to_year, values, analysis_dimension,
                              estimator = NA_character_, fixed_K = NA_integer_) {
  channels <- names(values)
  data.frame(
    from_year = as.integer(from_year),
    to_year = as.integer(to_year),
    transition = paste0(from_year, "-", to_year),
    analysis_dimension = analysis_dimension,
    estimator = estimator,
    fixed_K = fixed_K,
    family = vapply(channels, metric_family, character(1)),
    scale = vapply(channels, metric_scale, character(1)),
    channel = channels,
    score = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

for (i in seq_len(length(years) - 1L)) {
  from_year <- years[i]
  to_year <- years[i + 1L]
  A1 <- adjacency[[from_year]]
  A2 <- adjacency[[to_year]]
  before_indicators <- mx_network_indicators(A1)
  after_indicators <- mx_network_indicators(A2)
  transition_rows[[row_id]] <- append_transition(
    from_year, to_year,
    mx_indicator_change_scores(before_indicators, after_indicators),
    analysis_dimension = "conventional"
  )
  row_id <- row_id + 1L

  reference_pair <- mx_pair_scores(
    fits$fixed_estimator$block[[from_year]],
    fits$fixed_estimator$block[[to_year]],
    A1, A2,
    cut_samples = MX_CONFIG$empirical_cut_samples,
    cut_restarts = MX_CONFIG$cut_restarts
  )
  benchmark_values <- reference_pair[grepl("^singular_", names(reference_pair))]
  transition_rows[[row_id]] <- append_transition(
    from_year, to_year, benchmark_values,
    analysis_dimension = "adjacency benchmark"
  )
  row_id <- row_id + 1L

  for (method in MX_CONFIG$estimators) {
    g1 <- fits$fixed_estimator[[method]][[from_year]]
    g2 <- fits$fixed_estimator[[method]][[to_year]]
    values <- mx_pair_scores(
      g1, g2, A1, A2,
      cut_samples = MX_CONFIG$empirical_cut_samples,
      cut_restarts = MX_CONFIG$cut_restarts
    )
    values <- values[grepl("^graphon_", names(values))]
    values <- c(
      values,
      estimated_role_instability = 1 - mx_adjusted_rand_index(
        g1$node_labels, g2$node_labels
      )
    )
    transition_rows[[row_id]] <- append_transition(
      from_year, to_year, values,
      analysis_dimension = "estimator",
      estimator = method,
      fixed_K = MX_CONFIG$reference_K
    )
    row_id <- row_id + 1L
  }

  for (K in MX_CONFIG$resolution_grid) {
    key <- paste0("K", K)
    g1 <- fits$resolution[[key]][[from_year]]
    g2 <- fits$resolution[[key]][[to_year]]
    values <- mx_pair_scores(
      g1, g2, A1, A2,
      cut_samples = MX_CONFIG$empirical_cut_samples,
      cut_restarts = MX_CONFIG$cut_restarts
    )
    values <- values[grepl("^graphon_", names(values))]
    values <- c(
      values,
      estimated_role_instability = 1 - mx_adjusted_rand_index(
        g1$node_labels, g2$node_labels
      )
    )
    transition_rows[[row_id]] <- append_transition(
      from_year, to_year, values,
      analysis_dimension = "resolution",
      estimator = "block",
      fixed_K = K
    )
    row_id <- row_id + 1L
  }
}

transition_scores <- do.call(rbind, transition_rows)
group_estimator <- ifelse(
  is.na(transition_scores$estimator) | transition_scores$estimator == "",
  "not_applicable", transition_scores$estimator
)
group_K <- ifelse(
  is.na(transition_scores$fixed_K),
  "not_applicable", as.character(transition_scores$fixed_K)
)
transition_scores$evidence_percentile <- ave(
  transition_scores$score,
  interaction(
    transition_scores$analysis_dimension,
    group_estimator,
    group_K,
    transition_scores$channel,
    drop = TRUE,
    lex.order = TRUE
  ),
  FUN = mx_percent_rank
)
transition_scores$robust_z <- ave(
  transition_scores$score,
  interaction(
    transition_scores$analysis_dimension,
    group_estimator,
    group_K,
    transition_scores$channel,
    drop = TRUE,
    lex.order = TRUE
  ),
  FUN = function(x) {
    med <- stats::median(x, na.rm = TRUE)
    spread <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(spread) || spread <= 0) return(rep(0, length(x)))
    (x - med) / spread
  }
)
mx_write_csv(
  transition_scores,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_transition_evidence.csv")
)

agreement_table <- function(data, comparison_variable, reference_value = NULL) {
  channels <- sort(unique(data$channel))
  out <- list()
  k <- 1L
  for (channel in channels) {
    d <- data[data$channel == channel, , drop = FALSE]
    levels_available <- names(split(d$score, d[[comparison_variable]]))
    if (is.null(reference_value)) {
      pairs <- utils::combn(levels_available, 2L, simplify = FALSE)
    } else {
      reference_value <- as.character(reference_value)
      targets <- setdiff(levels_available, reference_value)
      pairs <- lapply(targets, function(x) c(reference_value, x))
    }
    for (pair in pairs) {
      d1 <- d[d[[comparison_variable]] == pair[1], c("transition", "score")]
      d2 <- d[d[[comparison_variable]] == pair[2], c("transition", "score")]
      m <- merge(d1, d2, by = "transition", suffixes = c("_1", "_2"))
      names1 <- setNames(m$score_1, m$transition)
      names2 <- setNames(m$score_2, m$transition)
      out[[k]] <- data.frame(
        channel = channel,
        comparison_1 = pair[1],
        comparison_2 = pair[2],
        spearman_rho = mx_spearman(m$score_1, m$score_2),
        top5_jaccard = mx_jaccard(mx_top_set(names1), mx_top_set(names2)),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

estimator_data <- transition_scores[
  transition_scores$analysis_dimension == "estimator", , drop = FALSE
]
estimator_agreement <- agreement_table(estimator_data, "estimator")
estimator_agreement$analysis <- "Estimator agreement at fixed K=13"
mx_write_csv(
  estimator_agreement,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_estimator_agreement.csv")
)

resolution_data <- transition_scores[
  transition_scores$analysis_dimension == "resolution", , drop = FALSE
]
resolution_data$fixed_K <- as.character(resolution_data$fixed_K)
resolution_agreement <- agreement_table(
  resolution_data, "fixed_K", reference_value = MX_CONFIG$reference_K
)
resolution_agreement$analysis <- "Resolution agreement versus K=13"
mx_write_csv(
  resolution_agreement,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_resolution_agreement.csv")
)

# Time-invariant destination reference set used only for interpretation.
incoming <- do.call(cbind, lapply(adjacency, colSums))
rownames(incoming) <- nodes$nodeID
core_definition <- nodes
core_definition$mean_in_strength_2006_2024 <- rowMeans(incoming)
core_order <- order(
  -core_definition$mean_in_strength_2006_2024,
  core_definition$nodeID
)
core_definition$destination_core_rank <- NA_integer_
core_definition$destination_core_rank[core_order] <- seq_len(nrow(core_definition))
core_definition$is_destination_core <-
  core_definition$destination_core_rank <= MX_CONFIG$destination_core_n
core_definition$is_regional_destination_core <-
  core_definition$is_destination_core & core_definition$nodeLabel != "Budapesti"
mx_write_csv(
  core_definition[order(core_definition$destination_core_rank), ],
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_destination_core_definition.csv")
)

role_memberships <- do.call(rbind, lapply(years, function(year) {
  mx_role_profile(
    fits$fixed_estimator$block[[year]],
    adjacency[[year]], nodes, year, core_definition
  )
}))
mx_write_csv(
  role_memberships,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_memberships_K13.csv")
)

role_year_summary <- do.call(rbind, lapply(years, function(year) {
  d <- role_memberships[role_memberships$year == as.integer(year), ]
  regional <- d[d$is_regional_destination_core, ]
  core <- d[d$is_destination_core, ]
  data.frame(
    year = as.integer(year),
    fixed_K = MX_CONFIG$reference_K,
    destination_core_effective_roles = mx_entropy_effective_n(table(core$destination_role_rank)),
    destination_core_same_role_rate = mx_same_group_rate(core$destination_role_rank),
    destination_core_pairwise_rank_gap = mx_mean_pairwise_abs(core$fitted_in_rank01),
    regional_core_effective_roles = mx_entropy_effective_n(table(regional$destination_role_rank)),
    regional_core_same_role_rate = mx_same_group_rate(regional$destination_role_rank),
    regional_core_pairwise_rank_gap = mx_mean_pairwise_abs(regional$fitted_in_rank01),
    stringsAsFactors = FALSE
  )
}))
mx_write_csv(
  role_year_summary,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_differentiation_by_year.csv")
)

role_transitions <- do.call(rbind, lapply(seq_len(length(years) - 1L), function(i) {
  y1 <- as.integer(years[i])
  y2 <- as.integer(years[i + 1L])
  a <- role_memberships[role_memberships$year == y1, ]
  b <- role_memberships[role_memberships$year == y2, ]
  m <- merge(
    a[, c("nodeID", "block", "fitted_in_rank01", "is_destination_core",
          "is_regional_destination_core")],
    b[, c("nodeID", "block", "fitted_in_rank01")],
    by = "nodeID", suffixes = c("_from", "_to"), sort = FALSE
  )
  core <- m$is_destination_core
  regional <- m$is_regional_destination_core
  data.frame(
    from_year = y1,
    to_year = y2,
    transition = paste0(y1, "-", y2),
    ARI_all_nodes = mx_adjusted_rand_index(m$block_from, m$block_to),
    ARI_destination_core = mx_adjusted_rand_index(
      m$block_from[core], m$block_to[core]
    ),
    ARI_regional_core = mx_adjusted_rand_index(
      m$block_from[regional], m$block_to[regional]
    ),
    mean_abs_fitted_in_rank_shift_all = mean(
      abs(m$fitted_in_rank01_to - m$fitted_in_rank01_from)
    ),
    mean_abs_fitted_in_rank_shift_regional = mean(
      abs(m$fitted_in_rank01_to[regional] - m$fitted_in_rank01_from[regional])
    ),
    stringsAsFactors = FALSE
  )
}))
mx_write_csv(
  role_transitions,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_transition_summary.csv")
)

refinement_rows <- list()
r <- 1L
for (year in as.character(MX_CONFIG$multiresolution_years)) {
  role_by_K <- lapply(MX_CONFIG$resolution_grid, function(K) {
    profile <- mx_role_profile(
      fits$resolution[[paste0("K", K)]][[year]],
      adjacency[[year]], nodes, year, core_definition
    )
    profile$destination_role_rank
  })
  names(role_by_K) <- paste0("K", MX_CONFIG$resolution_grid)
  for (j in seq_len(length(MX_CONFIG$resolution_grid) - 1L)) {
    low_K <- MX_CONFIG$resolution_grid[j]
    high_K <- MX_CONFIG$resolution_grid[j + 1L]
    refinement_rows[[r]] <- mx_refinement_summary(
      role_by_K[[paste0("K", low_K)]],
      role_by_K[[paste0("K", high_K)]],
      low_K, high_K, year
    )
    r <- r + 1L
  }
}
refinement_summary <- do.call(rbind, refinement_rows)
mx_write_csv(
  refinement_summary,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_multiresolution_refinement.csv")
)

flow_year <- as.character(MX_CONFIG$role_flow_year)
flow_roles <- role_memberships[
  role_memberships$year == MX_CONFIG$role_flow_year, , drop = FALSE
]
flow_roles <- flow_roles[match(nodes$nodeID, flow_roles$nodeID), ]
A_flow <- adjacency[[flow_year]]
flow_grid <- expand.grid(
  origin_role_rank = seq_len(MX_CONFIG$reference_K),
  destination_role_rank = seq_len(MX_CONFIG$reference_K),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
flow_grid$observed_weight <- vapply(seq_len(nrow(flow_grid)), function(i) {
  from <- flow_roles$destination_role_rank == flow_grid$origin_role_rank[i]
  to <- flow_roles$destination_role_rank == flow_grid$destination_role_rank[i]
  sum(A_flow[from, to, drop = FALSE])
}, numeric(1))
flow_grid$observed_weight_share <- flow_grid$observed_weight / sum(flow_grid$observed_weight)
flow_grid$year <- MX_CONFIG$role_flow_year
mx_write_csv(
  flow_grid,
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_flow_matrix_2018.csv")
)

saveRDS(
  list(
    transition_scores = transition_scores,
    estimator_agreement = estimator_agreement,
    resolution_agreement = resolution_agreement,
    core_definition = core_definition,
    role_memberships = role_memberships,
    role_year_summary = role_year_summary,
    role_transitions = role_transitions,
    refinement_summary = refinement_summary,
    role_flow = flow_grid
  ),
  file.path(MX_PATHS$object_output, "EMPIRICAL_analysis_results.rds"),
  compress = "gzip"
)

mx_log("Empirical benchmark completed.")



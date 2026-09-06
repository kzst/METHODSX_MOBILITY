# Known-truth simulation benchmark for calibrated transition screening.

mx_log("Stage 04: running the known-truth simulation benchmark.")

local({
  summary_path <- file.path(
    MX_PATHS$derived_simulation, "SIMULATION_detection_summary.csv"
  )
  object_path <- file.path(
    MX_PATHS$object_output, "SIMULATION_benchmark_results.rds"
  )

  if (MX_CONFIG$skip_simulation) {
    mx_assert_files(c(
      summary_path,
      file.path(
        MX_PATHS$derived_simulation, "SIMULATION_scenario_specification.csv"
      ),
      file.path(MX_PATHS$derived_simulation, "SIMULATION_kernel_design.csv")
    ), "Precomputed simulation output")
    mx_log("Simulation stage skipped; existing summary retained.")
    return(invisible(NULL))
  }

  signature <- list(
    seed = MX_CONFIG$seed,
    n = MX_CONFIG$simulation_n,
    K = MX_CONFIG$simulation_K,
    replicates = MX_CONFIG$simulation_replicates,
    cut_algorithm = MX_CONFIG$cut_algorithm,
    cut_samples = MX_CONFIG$simulation_cut_samples,
    cut_restarts = MX_CONFIG$cut_restarts,
    effects = MX_CONFIG$simulation_effects
  )
  if (MX_CONFIG$reuse_intermediate && file.exists(object_path)) {
    existing <- readRDS(object_path)
    if (identical(existing$signature, signature)) {
      reusable_outputs <- c(
        summary_path,
        file.path(
          MX_PATHS$derived_simulation, "SIMULATION_scenario_specification.csv"
        ),
        file.path(
          MX_PATHS$derived_simulation, "SIMULATION_kernel_design.csv"
        )
      )
      if (all(file.exists(reusable_outputs))) {
        mx_log("Reusing simulation benchmark from ", object_path, ".")
        return(invisible(NULL))
      }
    }
  }

  scenario_specs <- data.frame(
    scenario = c(
      "null", "global_scale", "density_shift", "destination_attraction",
      "directional_shift", "local_block_shock", "role_reassignment"
    ),
    scenario_label = c(
      "Null: sampling variation", "Global volume scaling",
      "Global density shift", "Destination-attraction shock",
      "Directional asymmetry shift", "Localized block-pair shock",
      "Node role reassignment"
    ),
    structural_scope = c(
      "none", "global magnitude only", "global topology",
      "one destination role", "selected directed block pairs",
      "one directed block pair", "subset of nodes"
    ),
    known_change = c(
      "No parameter change", "All conditional weight means increase",
      "All edge probabilities increase",
      "Incoming probability and weight increase for the leading destination block",
      "Selected directions increase while their reverse directions decrease",
      "One origin-to-destination block pair strengthens",
      "A controlled fraction of nodes changes latent role"
    ),
    diagnostic_purpose = c(
      "False-positive calibration", "Separate scale from structural shape",
      "Detect sparsity and connectivity change", "Detect target-specific attraction",
      "Detect directed reorganization", "Detect localized structural change",
      "Detect partition and role change"
    ),
    stringsAsFactors = FALSE
  )
  mx_write_csv(
    scenario_specs,
    file.path(MX_PATHS$derived_simulation, "SIMULATION_scenario_specification.csv")
  )

  K <- MX_CONFIG$simulation_K
  n <- MX_CONFIG$simulation_n
  base_labels <- rep(seq_len(K), length.out = n)
  base_probability <- matrix(c(
    0.50, 0.22, 0.18, 0.12,
    0.30, 0.45, 0.20, 0.16,
    0.24, 0.28, 0.40, 0.20,
    0.18, 0.22, 0.30, 0.35
  ), nrow = K, byrow = TRUE)
  base_mean <- matrix(c(
    5.0, 3.8, 2.4, 1.6,
    3.0, 4.5, 2.8, 1.9,
    2.2, 2.7, 4.0, 2.5,
    1.8, 2.1, 3.2, 3.6
  ), nrow = K, byrow = TRUE)

  scenario_parameters <- function(scenario, effect, labels, seed) {
    probability <- base_probability
    conditional_mean <- base_mean
    changed_labels <- labels
    if (scenario == "global_scale") {
      conditional_mean <- conditional_mean * (1 + 1.50 * effect)
    } else if (scenario == "density_shift") {
      probability <- probability * (1 + 0.90 * effect)
      probability[probability > 0.95] <- 0.95
    } else if (scenario == "destination_attraction") {
      probability[, 1L] <- pmin(0.95, probability[, 1L] * (1 + 0.60 * effect))
      conditional_mean[, 1L] <- conditional_mean[, 1L] * (1 + 1.50 * effect)
    } else if (scenario == "directional_shift") {
      increase <- rbind(c(1L, 2L), c(3L, 1L))
      decrease <- rbind(c(2L, 1L), c(1L, 3L))
      for (i in seq_len(nrow(increase))) {
        a <- increase[i, 1L]
        b <- increase[i, 2L]
        probability[a, b] <- min(0.95, probability[a, b] * (1 + effect))
        conditional_mean[a, b] <- conditional_mean[a, b] * (1 + 1.50 * effect)
        a <- decrease[i, 1L]
        b <- decrease[i, 2L]
        probability[a, b] <- max(0.01, probability[a, b] * (1 - 0.50 * effect))
        conditional_mean[a, b] <- max(
          0.05, conditional_mean[a, b] * (1 - 0.60 * effect)
        )
      }
    } else if (scenario == "local_block_shock") {
      probability[3L, 4L] <- min(
        0.95, probability[3L, 4L] * (1 + 1.50 * effect)
      )
      conditional_mean[3L, 4L] <-
        conditional_mean[3L, 4L] * (1 + 2.00 * effect)
    } else if (scenario == "role_reassignment") {
      candidates <- which(labels == K)
      fraction <- min(0.60, 0.10 + 0.40 * effect)
      count <- max(1L, round(length(candidates) * fraction))
      set.seed(seed)
      moved <- sample(candidates, count, replace = FALSE)
      changed_labels[moved] <- 2L
    } else if (scenario != "null") {
      stop("Unknown simulation scenario: ", scenario)
    }
    if (
      !is.matrix(probability) ||
        !identical(dim(probability), c(K, K)) ||
        any(!is.finite(probability)) ||
        any(probability < 0 | probability > 1)
    ) {
      stop("Scenario probability parameters must form a finite K x K matrix in [0, 1].")
    }
    if (
      !is.matrix(conditional_mean) ||
        !identical(dim(conditional_mean), c(K, K)) ||
        any(!is.finite(conditional_mean)) ||
        any(conditional_mean < 0)
    ) {
      stop("Scenario weight parameters must form a finite nonnegative K x K matrix.")
    }
    list(
      probability = probability,
      conditional_mean = conditional_mean,
      labels = changed_labels
    )
  }

  simulate_network <- function(parameters, theta_out, theta_in, seed) {
    labels <- parameters$labels
    probability <- parameters$probability[labels, labels, drop = FALSE]
    lambda <- parameters$conditional_mean[labels, labels, drop = FALSE] *
      outer(theta_out, theta_in)
    set.seed(seed)
    present <- matrix(
      stats::rbinom(n * n, 1L, as.numeric(probability)),
      nrow = n, ncol = n
    )
    weights <- matrix(
      stats::rpois(n * n, lambda = as.numeric(lambda)),
      nrow = n, ncol = n
    )
    A <- present * (1 + weights)
    diag(A) <- 0
    storage.mode(A) <- "double"
    A
  }

  expected_network <- function(parameters) {
    labels <- parameters$labels
    probability <- parameters$probability[labels, labels, drop = FALSE]
    conditional_mean <- parameters$conditional_mean[labels, labels, drop = FALSE]
    E <- probability * (1 + conditional_mean)
    diag(E) <- 0
    E
  }

  aggregate_by_reference_blocks <- function(M, labels) {
    out <- matrix(0, K, K)
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        out[i, j] <- mean(M[labels == i, labels == j, drop = FALSE])
      }
    }
    out
  }

  sim_family <- function(channel) {
    if (grepl("^indicator_", channel)) return("conventional indicator")
    if (grepl("^graphon_", channel)) return("fitted-kernel distance")
    if (grepl("^singular_", channel)) return("adjacency benchmark")
    if (channel == "estimated_role_instability") return("structural roles")
    "other"
  }

  sim_scale <- function(channel) {
    if (grepl("_shape$", channel)) return("unit-mean shape")
    if (grepl("_raw$", channel)) return("raw weighted scale")
    if (channel == "estimated_role_instability") return("label-invariant partition")
    "native indicator scale"
  }

  run_pair <- function(scenario, effect_name, effect, replicate_id) {
    scenario_index <- match(scenario, scenario_specs$scenario)
    effect_index <- match(effect_name, c("null", names(MX_CONFIG$simulation_effects)))
    pair_seed <- MX_CONFIG$seed + 100000L * scenario_index +
      1000L * effect_index + replicate_id
    set.seed(pair_seed)
    theta_out <- exp(stats::rnorm(n, mean = 0, sd = 0.30))
    theta_in <- exp(stats::rnorm(n, mean = 0, sd = 0.30))
    theta_out <- theta_out / mean(theta_out)
    theta_in <- theta_in / mean(theta_in)

    baseline <- scenario_parameters("null", 0, base_labels, pair_seed)
    changed <- scenario_parameters(scenario, effect, base_labels, pair_seed + 1L)
    A0 <- simulate_network(baseline, theta_out, theta_in, pair_seed + 2L)
    A1 <- simulate_network(changed, theta_out, theta_in, pair_seed + 3L)
    g0 <- mx_fit_graphon(A0, K = K, method = "block")
    g1 <- mx_fit_graphon(A1, K = K, method = "block")

    indicator_scores <- mx_indicator_change_scores(
      mx_network_indicators(A0), mx_network_indicators(A1)
    )
    structural_scores <- mx_pair_scores(
      g0, g1, A0, A1,
      cut_samples = MX_CONFIG$simulation_cut_samples,
      cut_restarts = MX_CONFIG$cut_restarts
    )
    scores <- c(
      indicator_scores,
      structural_scores,
      estimated_role_instability = 1 - mx_adjusted_rand_index(
        g0$node_labels, g1$node_labels
      )
    )
    data.frame(
      scenario = scenario,
      scenario_label = scenario_specs$scenario_label[scenario_index],
      effect_label = effect_name,
      effect = effect,
      replicate = replicate_id,
      seed = pair_seed,
      known_role_instability = 1 - mx_adjusted_rand_index(
        baseline$labels, changed$labels
      ),
      family = vapply(names(scores), sim_family, character(1)),
      scale = vapply(names(scores), sim_scale, character(1)),
      channel = names(scores),
      score = as.numeric(scores),
      stringsAsFactors = FALSE
    )
  }

  run_design <- data.frame(
    scenario = "null",
    effect_label = "null",
    effect = 0,
    stringsAsFactors = FALSE
  )
  for (scenario in setdiff(scenario_specs$scenario, "null")) {
    run_design <- rbind(
      run_design,
      data.frame(
        scenario = scenario,
        effect_label = names(MX_CONFIG$simulation_effects),
        effect = as.numeric(MX_CONFIG$simulation_effects),
        stringsAsFactors = FALSE
      )
    )
  }

  raw_rows <- vector(
    "list", nrow(run_design) * MX_CONFIG$simulation_replicates
  )
  row_id <- 1L
  for (design_id in seq_len(nrow(run_design))) {
    design_row <- run_design[design_id, ]
    mx_log(
      "Simulation: ", design_row$scenario, ", effect=", design_row$effect_label,
      ", replicates=", MX_CONFIG$simulation_replicates, "."
    )
    for (replicate_id in seq_len(MX_CONFIG$simulation_replicates)) {
      raw_rows[[row_id]] <- run_pair(
        design_row$scenario,
        design_row$effect_label,
        design_row$effect,
        replicate_id
      )
      row_id <- row_id + 1L
    }
  }
  raw_results <- do.call(rbind, raw_rows)

  null_results <- raw_results[raw_results$scenario == "null", ]
  threshold_rows <- lapply(split(null_results, null_results$channel), function(d) {
    data.frame(
      channel = d$channel[1L],
      null_threshold_95 = as.numeric(
        stats::quantile(d$score, probs = 0.95, na.rm = TRUE, names = FALSE)
      ),
      null_median = stats::median(d$score, na.rm = TRUE),
      null_mad = stats::mad(d$score, constant = 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  thresholds <- do.call(rbind, threshold_rows)
  raw_results <- merge(raw_results, thresholds, by = "channel", all.x = TRUE)
  raw_results$detected <- raw_results$score > raw_results$null_threshold_95

  summary_groups <- split(
    raw_results,
    interaction(
      raw_results$scenario,
      raw_results$effect_label,
      raw_results$channel,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  detection_summary <- do.call(rbind, lapply(summary_groups, function(d) {
    data.frame(
      scenario = d$scenario[1L],
      scenario_label = d$scenario_label[1L],
      effect_label = d$effect_label[1L],
      effect = d$effect[1L],
      family = d$family[1L],
      scale = d$scale[1L],
      channel = d$channel[1L],
      replicates = nrow(d),
      calibrated_detection_rate = mean(d$detected, na.rm = TRUE),
      mean_score = mean(d$score, na.rm = TRUE),
      median_score = stats::median(d$score, na.rm = TRUE),
      score_q025 = as.numeric(stats::quantile(d$score, 0.025, na.rm = TRUE)),
      score_q975 = as.numeric(stats::quantile(d$score, 0.975, na.rm = TRUE)),
      null_threshold_95 = d$null_threshold_95[1L],
      stringsAsFactors = FALSE
    )
  }))
  rownames(detection_summary) <- NULL

  kernel_rows <- list()
  krow <- 1L
  medium_effect <- unname(MX_CONFIG$simulation_effects["medium"])
  baseline_expected <- expected_network(
    scenario_parameters("null", 0, base_labels, MX_CONFIG$seed)
  )
  baseline_block <- aggregate_by_reference_blocks(baseline_expected, base_labels)
  for (scenario in setdiff(scenario_specs$scenario, "null")) {
    changed_parameters <- scenario_parameters(
      scenario, medium_effect, base_labels, MX_CONFIG$seed + match(
        scenario, scenario_specs$scenario
      )
    )
    changed_block <- aggregate_by_reference_blocks(
      expected_network(changed_parameters), base_labels
    )
    label <- scenario_specs$scenario_label[scenario_specs$scenario == scenario]
    for (state in c("Baseline", "Changed")) {
      M <- if (state == "Baseline") baseline_block else changed_block
      grid <- expand.grid(
        origin_block = seq_len(K),
        destination_block = seq_len(K),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      grid$expected_weight <- M[cbind(grid$origin_block, grid$destination_block)]
      grid$scenario <- scenario
      grid$scenario_label <- label
      grid$state <- state
      grid$effect_label <- "medium"
      grid$effect <- medium_effect
      kernel_rows[[krow]] <- grid
      krow <- krow + 1L
    }
  }
  kernel_design <- do.call(rbind, kernel_rows)

  mx_write_csv(
    raw_results,
    file.path(MX_PATHS$derived_simulation, "SIMULATION_raw_scores.csv")
  )
  mx_write_csv(
    thresholds,
    file.path(MX_PATHS$derived_simulation, "SIMULATION_null_thresholds.csv")
  )
  mx_write_csv(detection_summary, summary_path)
  mx_write_csv(
    kernel_design,
    file.path(MX_PATHS$derived_simulation, "SIMULATION_kernel_design.csv")
  )

  saveRDS(
    list(
      signature = signature,
      scenario_specs = scenario_specs,
      thresholds = thresholds,
      detection_summary = detection_summary,
      kernel_design = kernel_design
    ),
    object_path,
    compress = "gzip"
  )

  mx_log("Simulation benchmark completed.")
})



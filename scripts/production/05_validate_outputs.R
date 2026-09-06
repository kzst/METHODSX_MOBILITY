# Deterministic validation and reproducibility metadata.

mx_log("Stage 05: validating analytical outputs.")

required_outputs <- c(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_node_aligned_networks.rds"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_input_checks.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_indicators.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_transition_evidence.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_estimator_agreement.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_resolution_agreement.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_memberships_K13.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_transition_summary.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_multiresolution_refinement.csv"),
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_flow_matrix_2018.csv"),
  file.path(MX_PATHS$derived_simulation, "SIMULATION_scenario_specification.csv"),
  file.path(MX_PATHS$derived_simulation, "SIMULATION_detection_summary.csv"),
  file.path(MX_PATHS$derived_simulation, "SIMULATION_kernel_design.csv")
)
mx_assert_files(required_outputs, "Pipeline output")

input_checks <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_input_checks.csv")
)
annual_indicators <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_annual_indicators.csv")
)
transition_evidence <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_transition_evidence.csv")
)
role_memberships <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_memberships_K13.csv")
)
refinement <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_multiresolution_refinement.csv")
)
role_flow <- utils::read.csv(
  file.path(MX_PATHS$derived_empirical, "EMPIRICAL_role_flow_matrix_2018.csv")
)
simulation_summary <- utils::read.csv(
  file.path(MX_PATHS$derived_simulation, "SIMULATION_detection_summary.csv")
)

expected_indicators <- c(
  "total_weight", "edge_density", "weighted_reciprocity",
  "in_strength_gini", "out_strength_gini", "effective_destinations",
  "top10_destination_share", "binary_transitivity"
)
expected_scenarios <- c(
  "null", "global_scale", "density_shift", "destination_attraction",
  "directional_shift", "local_block_shock", "role_reassignment"
)
simulation_rates_valid <-
  all(is.finite(simulation_summary$calibrated_detection_rate)) &&
  all(
    simulation_summary$calibrated_detection_rate >= 0 &
      simulation_summary$calibrated_detection_rate <= 1
  )
null_rates <- simulation_summary$calibrated_detection_rate[
  simulation_summary$scenario == "null"
]
maximum_null_rate <- if (length(null_rates) && all(is.finite(null_rates))) {
  max(null_rates)
} else {
  NA_real_
}

checks <- data.frame(
  check = c(
    "pipeline_version",
    "years_2006_2024",
    "nodes_175_each_year",
    "year_2016_total_weight",
    "year_2019_total_weight",
    "annual_indicator_coverage",
    "eighteen_consecutive_transitions",
    "finite_transition_scores",
    "three_estimators_present",
    "four_resolutions_present",
    "role_membership_dimensions",
    "multiresolution_refinement_dimensions",
    "role_flow_matrix_dimensions",
    "simulation_scenarios_present",
    "simulation_detection_rates_valid",
    "null_false_positive_rate_calibrated"
  ),
  observed = c(
    MX_CONFIG$pipeline_version,
    paste(input_checks$year, collapse = ","),
    paste(sort(unique(input_checks$nodes)), collapse = ","),
    as.character(input_checks$total_weight[input_checks$year == 2016]),
    as.character(input_checks$total_weight[input_checks$year == 2019]),
    as.character(nrow(annual_indicators)),
    as.character(length(unique(transition_evidence$transition))),
    as.character(all(is.finite(transition_evidence$score))),
    paste(sort(unique(transition_evidence$estimator[
      transition_evidence$analysis_dimension == "estimator"
    ])), collapse = ","),
    paste(sort(unique(transition_evidence$fixed_K[
      transition_evidence$analysis_dimension == "resolution"
    ])), collapse = ","),
    paste(nrow(role_memberships), ncol(role_memberships), sep = "x"),
    as.character(nrow(refinement)),
    paste(nrow(role_flow), ncol(role_flow), sep = "x"),
    paste(sort(unique(simulation_summary$scenario)), collapse = ","),
    as.character(simulation_rates_valid),
    as.character(maximum_null_rate)
  ),
  expected = c(
    "0.2.0",
    paste(2006:2024, collapse = ","),
    "175",
    "358367",
    "384035",
    as.character(19L * length(expected_indicators)),
    "18",
    "TRUE",
    paste(sort(MX_CONFIG$estimators), collapse = ","),
    paste(sort(MX_CONFIG$resolution_grid), collapse = ","),
    paste(19L * MX_CONFIG$expected_nodes, 18L, sep = "x"),
    as.character(
      length(MX_CONFIG$multiresolution_years) *
        (length(MX_CONFIG$resolution_grid) - 1L)
    ),
    paste(MX_CONFIG$reference_K^2, 5L, sep = "x"),
    paste(sort(expected_scenarios), collapse = ","),
    "TRUE",
    "<=0.10"
  ),
  stringsAsFactors = FALSE
)

checks$passed <- c(
  identical(MX_CONFIG$pipeline_version, "0.2.0"),
  identical(as.integer(input_checks$year), as.integer(2006:2024)),
  all(input_checks$nodes == MX_CONFIG$expected_nodes),
  isTRUE(all.equal(
    as.numeric(input_checks$total_weight[input_checks$year == 2016]),
    358367,
    tolerance = 0
  )),
  isTRUE(all.equal(
    as.numeric(input_checks$total_weight[input_checks$year == 2019]),
    384035,
    tolerance = 0
  )),
  nrow(annual_indicators) == 19L * length(expected_indicators) &&
    identical(sort(unique(annual_indicators$indicator)), sort(expected_indicators)),
  length(unique(transition_evidence$transition)) == 18L,
  all(is.finite(transition_evidence$score)),
  identical(
    sort(unique(transition_evidence$estimator[
      transition_evidence$analysis_dimension == "estimator"
    ])),
    sort(MX_CONFIG$estimators)
  ),
  identical(
    sort(unique(as.integer(transition_evidence$fixed_K[
      transition_evidence$analysis_dimension == "resolution"
    ]))),
    sort(MX_CONFIG$resolution_grid)
  ),
  nrow(role_memberships) == 19L * MX_CONFIG$expected_nodes &&
    ncol(role_memberships) == 18L,
  nrow(refinement) == length(MX_CONFIG$multiresolution_years) *
    (length(MX_CONFIG$resolution_grid) - 1L),
  nrow(role_flow) == MX_CONFIG$reference_K^2 &&
    ncol(role_flow) == 5L &&
    isTRUE(all.equal(sum(role_flow$observed_weight_share), 1, tolerance = 1e-12)),
  identical(
    sort(unique(simulation_summary$scenario)),
    sort(expected_scenarios)
  ),
  simulation_rates_valid,
  is.finite(maximum_null_rate) && maximum_null_rate <= 0.10
)

validation_path <- file.path(
  MX_PATHS$diagnostics, "METHODSX_PIPELINE_validation.csv"
)
mx_write_csv(checks, validation_path)

config_table <- data.frame(
  setting = c(
    "pipeline_version", "mode", "seed", "reference_K", "resolution_grid",
    "estimators", "cut_algorithm", "cut_restarts",
    "empirical_cut_samples", "simulation_replicates",
    "simulation_cut_samples", "simulation_n", "simulation_K"
  ),
  value = c(
    MX_CONFIG$pipeline_version,
    MX_CONFIG$mode,
    MX_CONFIG$seed,
    MX_CONFIG$reference_K,
    paste(MX_CONFIG$resolution_grid, collapse = ";"),
    paste(MX_CONFIG$estimators, collapse = ";"),
    MX_CONFIG$cut_algorithm,
    MX_CONFIG$cut_restarts,
    MX_CONFIG$empirical_cut_samples,
    MX_CONFIG$simulation_replicates,
    MX_CONFIG$simulation_cut_samples,
    MX_CONFIG$simulation_n,
    MX_CONFIG$simulation_K
  ),
  stringsAsFactors = FALSE
)
mx_write_csv(
  config_table,
  file.path(MX_PATHS$diagnostics, "METHODSX_PIPELINE_configuration.csv")
)

mx_refresh_reproducibility_metadata <- function(include_rendered_outputs = FALSE) {
  session_path <- file.path(
    MX_PATHS$diagnostics, "METHODSX_PIPELINE_sessionInfo.txt"
  )
  writeLines(capture.output(sessionInfo()), session_path)

  diagnostic_files <- file.path(
    MX_PATHS$diagnostics,
    c(
      "METHODSX_PIPELINE_validation.csv",
      "METHODSX_PIPELINE_configuration.csv",
      "METHODSX_PIPELINE_sessionInfo.txt"
    )
  )
  generated <- unique(c(
    list.files("data/derived", recursive = TRUE, full.names = TRUE),
    list.files(MX_PATHS$object_output, recursive = TRUE, full.names = TRUE),
    diagnostic_files
  ))

  if (isTRUE(include_rendered_outputs)) {
    generated <- unique(c(
      generated,
      list.files(MX_PATHS$figure_output, recursive = TRUE, full.names = TRUE),
      file.path(MX_PATHS$diagnostics, "METHODSX_FIGURE_manifest.csv"),
      file.path("output", c("MethodsX_MOB_figures.html", "MethodsX_MOB.docx"))
    ))
  }

  manifest_path <- file.path(
    MX_PATHS$diagnostics, "METHODSX_PIPELINE_output_manifest.csv"
  )
  generated <- sort(unique(gsub("\\\\", "/", generated)))
  generated <- setdiff(generated, gsub("\\\\", "/", manifest_path))
  generated <- generated[file.exists(generated) & !dir.exists(generated)]

  manifest <- data.frame(
    file = generated,
    bytes = unname(file.info(generated)$size),
    modified = as.character(file.info(generated)$mtime),
    md5 = unname(tools::md5sum(generated)),
    stringsAsFactors = FALSE
  )
  mx_write_csv(manifest, manifest_path)
  invisible(manifest_path)
}

print(checks, row.names = FALSE)
if (!all(checks$passed)) {
  stop(
    "MethodsX pipeline validation failed. Review: ", validation_path,
    call. = FALSE
  )
}

mx_refresh_reproducibility_metadata(include_rendered_outputs = FALSE)
mx_log("All analytical output checks passed.")




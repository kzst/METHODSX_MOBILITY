# Standalone release-candidate checks.

required_sources <- c(
  "METHODSX_MOB.Rproj",
  "MethodsX_MOB.Rmd",
  "MethodsX-Method-Article-Template.docx",
  "run_all.R",
  "DESCRIPTION",
  "CITATION.cff",
    file.path("scripts", "production", sprintf("%02d_%s.R", 0:9, c(
    "config", "prepare_empirical", "analysis_helpers", "empirical_benchmark",
    "simulation_benchmark", "validate_outputs", "graphical_abstract",
    "freeze_docx_fields", "cover_letter", "validate_template_compliance"
  )))
)

missing_sources <- required_sources[!file.exists(required_sources)]
if (length(missing_sources)) {
  stop(
    "Release source file(s) missing: ",
    paste(missing_sources, collapse = ", "),
    call. = FALSE
  )
}

edge_files <- list.files(
  file.path("data", "raw", "empirical"),
  pattern = "^[0-9]{4}_from_to_kist_weight[.]txt$",
  full.names = TRUE
)
edge_years <- sort(suppressWarnings(as.integer(substr(basename(edge_files), 1, 4))))
if (!identical(edge_years, 2006:2024)) {
  stop("Raw empirical files do not cover 2006–2024 exactly.", call. = FALSE)
}

validation_path <- file.path(
  "output", "diagnostics", "METHODSX_PIPELINE_validation.csv"
)
if (!file.exists(validation_path)) {
  stop("Validation report missing; run source('run_all.R') first.", call. = FALSE)
}
validation <- utils::read.csv(validation_path, stringsAsFactors = FALSE)
if (nrow(validation) != 16L || !all(validation$passed)) {
  stop("The 16 analytical validation checks have not all passed.", call. = FALSE)
}

require_rendered <- tolower(Sys.getenv(
  "METHODSX_REQUIRE_RENDERED", unset = "1"
)) %in% c("1", "true", "yes", "y", "on")

if (require_rendered) {
  figure_stems <- c(
    "FIG00_graphical_abstract", "FIG01_workflow",
    "FIG02_simulation_scenarios", "FIG03_simulation_power",
    "FIG04_empirical_transition_evidence", "FIG05_robustness_agreement",
    "FIG06_role_flow_matrix", "FIGS01_refinement_consistency"
  )
  expected_figures <- as.vector(outer(
    file.path("output", "figures", figure_stems),
    c(".pdf", ".png"), paste0
  ))
  required_rendered <- c(
    expected_figures,
    file.path("output", "MethodsX_MOB_figures.html"),
    file.path("output", "MethodsX_MOB.docx"),
    file.path("output", "Cover_Letter_MethodsX.docx"),
    file.path("submission", "Cover_Letter_MethodsX.docx"),
    file.path("output", "diagnostics", "METHODSX_FIGURE_manifest.csv"),
    file.path("output", "diagnostics", "METHODSX_TEMPLATE_compliance.csv"),
    file.path("output", "diagnostics", "METHODSX_PIPELINE_output_manifest.csv")
  )
  missing_rendered <- required_rendered[!file.exists(required_rendered)]
  if (length(missing_rendered)) {
    stop(
      "Rendered release file(s) missing: ",
      paste(missing_rendered, collapse = ", "),
      call. = FALSE
    )
  }
  figure_manifest <- utils::read.csv(
    file.path("output", "diagnostics", "METHODSX_FIGURE_manifest.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(figure_manifest) != 16L) {
    stop("Figure manifest must contain 16 PDF/PNG exports.", call. = FALSE)
  }

  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required for the DOCX field audit.", call. = FALSE)
  }
  source(file.path("scripts", "production", "07_freeze_docx_fields.R"))
  docx_audit <- mx_audit_docx_fields(file.path("output", "MethodsX_MOB.docx"))
  source(file.path("scripts", "production", "09_validate_template_compliance.R"))
  mx_validate_methodsx_template(file.path("output", "MethodsX_MOB.docx"))
  if (
    docx_audit$seq_fields != 0L || docx_audit$field_markers != 0L ||
      docx_audit$dirty_fields != 0L || docx_audit$update_fields != 0L ||
      docx_audit$comment_anchors != 0L ||
      length(docx_audit$comment_parts) != 0L ||
      length(docx_audit$external_file_targets) != 0L ||
      length(docx_audit$invalid_relationship_ids) != 0L
  ) {
    stop(
      "The rendered Word manuscript still contains updateable fields or ",
      "external file targets, or the Word package contains invalid ",
      "relationship IDs.",
      call. = FALSE
    )
  }
}

message("Release checks passed.")

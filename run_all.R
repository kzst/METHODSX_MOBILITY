# Complete MethodsX MOB analytical pipeline.
#
# Diagnostic run (recommended first):
#   Sys.setenv(METHODSX_MODE = "diagnostic")
#   source("run_all.R")
#
# Publication run:
#   Sys.setenv(
#     METHODSX_MODE = "publication",
#     METHODSX_REUSE_INTERMEDIATE = "0"
#   )
#   source("run_all.R")

mx_find_project_root <- function() {
  if (file.exists("METHODSX_MOB.Rproj") && dir.exists("scripts/production")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    candidate <- dirname(normalizePath(
      sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = TRUE
    ))
    if (file.exists(file.path(candidate, "METHODSX_MOB.Rproj"))) return(candidate)
  }
  stop(
    "Run run_all.R from the METHODSX_MOB project root or open METHODSX_MOB.Rproj first.",
    call. = FALSE
  )
}

project_root <- mx_find_project_root()
setwd(project_root)

source(file.path("scripts", "production", "00_config.R"))
pipeline_started <- Sys.time()
mx_log("MethodsX MOB pipeline started in ", project_root, ".")

stages <- c(
  file.path("scripts", "production", "01_prepare_empirical.R"),
  file.path("scripts", "production", "02_analysis_helpers.R"),
  file.path("scripts", "production", "03_empirical_benchmark.R"),
  file.path("scripts", "production", "04_simulation_benchmark.R"),
  file.path("scripts", "production", "05_validate_outputs.R")
)
mx_assert_files(stages, "Production script")

for (stage in stages) {
  mx_log("Sourcing ", stage, ".")
  source(stage, local = .GlobalEnv, chdir = FALSE)
}

if (MX_CONFIG$render_rmd) {
  mx_require_packages(c(
    "rmarkdown", "knitr", "ggplot2", "dplyr", "tidyr", "scales",
    "officedown", "officer", "flextable", "xml2", "zip"
  ))
  mx_assert_files(
    c(
      "MethodsX_MOB.Rmd", "MethodsX-Method-Article-Template.docx",
      file.path("scripts", "production", c(
        "06_graphical_abstract.R", "07_freeze_docx_fields.R",
        "08_cover_letter.R", "09_validate_template_compliance.R"
      ))
    ),
    "Manuscript source"
  )
  source(
    file.path("scripts", "production", "07_freeze_docx_fields.R"),
    local = .GlobalEnv, chdir = FALSE
  )
  source(
    file.path("scripts", "production", "08_cover_letter.R"),
    local = .GlobalEnv, chdir = FALSE
  )
  source(
    file.path("scripts", "production", "09_validate_template_compliance.R"),
    local = .GlobalEnv, chdir = FALSE
  )

  mx_log("Rendering the MethodsX HTML review manuscript.")
  rmarkdown::render(
    input = "MethodsX_MOB.Rmd",
    output_format = "html_document",
    output_file = "MethodsX_MOB_figures.html",
    output_dir = "output",
    params = list(build_figures = TRUE),
    clean = TRUE,
    envir = new.env(parent = globalenv())
  )

  mx_log("Rendering the MethodsX Word manuscript.")
  rmarkdown::render(
    input = "MethodsX_MOB.Rmd",
    output_format = "officedown::rdocx_document",
    output_file = "MethodsX_MOB.docx",
    output_dir = "output",
    params = list(build_figures = TRUE),
    clean = TRUE,
    envir = new.env(parent = globalenv())
  )

  frozen_fields <- mx_freeze_docx_fields(
    file.path("output", "MethodsX_MOB.docx")
  )
  mx_log(
    "Converted ", frozen_fields,
    " automatic figure-number fields to fixed text in the Word manuscript."
  )
  compliance <- mx_validate_methodsx_template(
    file.path("output", "MethodsX_MOB.docx")
  )
  mx_log(
    "MethodsX template compliance passed: title ", compliance$title_words,
    " words; abstract ", compliance$abstract_words,
    " words; all required sections present in template order."
  )

  mx_log("Rendering the MethodsX cover letter.")
  mx_build_cover_letter(file.path("output", "Cover_Letter_MethodsX.docx"))
  dir.create("submission", recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(
    file.path("output", "Cover_Letter_MethodsX.docx"),
    file.path("submission", "Cover_Letter_MethodsX.docx"),
    overwrite = TRUE
  )) {
    stop("Could not refresh submission/Cover_Letter_MethodsX.docx.", call. = FALSE)
  }

  mx_assert_files(
    c(file.path("output", c(
      "MethodsX_MOB_figures.html", "MethodsX_MOB.docx",
      "Cover_Letter_MethodsX.docx"
    )), file.path(
      "output", "diagnostics", "METHODSX_TEMPLATE_compliance.csv"
    )),
    "Rendered manuscript"
  )

  figure_stems <- c(
    "FIG00_graphical_abstract",
    "FIG01_workflow",
    "FIG02_simulation_scenarios",
    "FIG03_simulation_power",
    "FIG04_empirical_transition_evidence",
    "FIG05_robustness_agreement",
    "FIG06_role_flow_matrix",
    "FIGS01_refinement_consistency"
  )
  expected_figures <- as.vector(outer(
    file.path(MX_PATHS$figure_output, figure_stems),
    c(".pdf", ".png"),
    paste0
  ))
  mx_assert_files(expected_figures, "Rendered figure")

  figure_manifest <- data.frame(
    file = gsub("\\\\", "/", expected_figures),
    bytes = unname(file.info(expected_figures)$size),
    md5 = unname(tools::md5sum(expected_figures)),
    stringsAsFactors = FALSE
  )
  mx_write_csv(
    figure_manifest,
    file.path(MX_PATHS$diagnostics, "METHODSX_FIGURE_manifest.csv")
  )
  mx_refresh_reproducibility_metadata(include_rendered_outputs = TRUE)
  mx_log(
    "All ", length(expected_figures),
    " PDF/PNG figure exports were verified."
  )
  mx_log("The final output manifest now includes the rendered manuscripts and figures.")
}

elapsed <- as.numeric(difftime(Sys.time(), pipeline_started, units = "mins"))
mx_log(
  "MethodsX MOB pipeline completed in ", sprintf("%.1f", elapsed), " minutes."
)
mx_log(
  "Validation report: output/diagnostics/METHODSX_PIPELINE_validation.csv"
)
if (MX_CONFIG$render_rmd) {
  mx_log("Figure review file: output/MethodsX_MOB_figures.html")
  mx_log("Template-compliant Word manuscript: output/MethodsX_MOB.docx")
  mx_log("Cover letter: output/Cover_Letter_MethodsX.docx")
}

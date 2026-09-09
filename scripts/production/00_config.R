# MethodsX MOB analysis configuration
# Version: 0.3.0

options(stringsAsFactors = FALSE)

mx_log <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

mx_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "1" else "0")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y", "on")
}

mx_env_integer <- function(name, default, minimum = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be one integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

mx_assert_files <- function(paths, context = "Required file") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(context, " missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(paths)
}

mx_require_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Install them with install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "), ")).",
      call. = FALSE
    )
  }
  invisible(packages)
}

mx_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create directory: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

mx_write_csv <- function(x, path) {
  mx_dir(dirname(path))
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

mx_percent_rank <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) == 1L) {
    out[ok] <- 1
  } else if (sum(ok) > 1L) {
    out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  }
  out
}

mx_mode <- tolower(trimws(Sys.getenv("METHODSX_MODE", unset = "diagnostic")))
if (!mx_mode %in% c("diagnostic", "publication")) {
  stop("METHODSX_MODE must be 'diagnostic' or 'publication'.", call. = FALSE)
}

mx_default_sim_reps <- if (mx_mode == "publication") 200L else 20L
mx_default_sim_cut <- if (mx_mode == "publication") 128L else 32L
mx_default_emp_cut <- if (mx_mode == "publication") 1000L else 100L
mx_default_cut_restarts <- if (mx_mode == "publication") 200L else 20L

MX_CONFIG <- list(
  pipeline_version = "0.3.0",
  mode = mx_mode,
  seed = 20260817L,
  years = 2006:2024,
  expected_nodes = 175L,
  reference_K = 13L,
  resolution_grid = c(2L, 4L, 8L, 13L),
  estimators = c("block", "smooth", "spectral"),
  destination_core_n = 20L,
  multiresolution_years = c(2015L, 2018L, 2021L, 2024L),
  role_flow_year = 2018L,
  cut_algorithm = "alternating_maximisation_lower_bound_v1",
  cut_restarts = mx_env_integer(
    "METHODSX_CUT_RESTARTS", mx_default_cut_restarts, minimum = 1L
  ),
  empirical_cut_samples = mx_env_integer(
    "METHODSX_EMPIRICAL_CUT_SAMPLES", mx_default_emp_cut, minimum = 10L
  ),
  simulation_replicates = mx_env_integer(
    "METHODSX_SIM_REPS", mx_default_sim_reps, minimum = 10L
  ),
  simulation_cut_samples = mx_env_integer(
    "METHODSX_SIM_CUT_SAMPLES", mx_default_sim_cut, minimum = 10L
  ),
  simulation_n = mx_env_integer("METHODSX_SIM_N", 80L, minimum = 40L),
  simulation_K = 4L,
  simulation_effects = c(small = 0.25, medium = 0.50, large = 1.00),
  reuse_intermediate = mx_env_flag("METHODSX_REUSE_INTERMEDIATE", default = TRUE),
  render_rmd = mx_env_flag("METHODSX_RENDER_RMD", default = TRUE),
  skip_simulation = mx_env_flag("METHODSX_SKIP_SIMULATION", default = FALSE)
)

MX_PATHS <- list(
  raw_empirical = file.path("data", "raw", "empirical"),
  derived_empirical = file.path("data", "derived", "empirical"),
  derived_simulation = file.path("data", "derived", "simulation"),
  figure_output = file.path("output", "figures"),
  object_output = file.path("output", "objects"),
  diagnostics = file.path("output", "diagnostics"),
  graphon_reference = file.path(
    "scripts", "reference", "revised", "graphon_distance_directed.R"
  )
)

invisible(lapply(MX_PATHS[c(
  "derived_empirical", "derived_simulation", "figure_output",
  "object_output", "diagnostics"
)], mx_dir))

mx_log(
  "Configuration loaded: mode=", MX_CONFIG$mode,
  ", simulation replicates=", MX_CONFIG$simulation_replicates,
  ", empirical cut-floor samples=", MX_CONFIG$empirical_cut_samples,
  ", simulation cut-floor samples=", MX_CONFIG$simulation_cut_samples,
  ", alternating cut restarts=", MX_CONFIG$cut_restarts, "."
)


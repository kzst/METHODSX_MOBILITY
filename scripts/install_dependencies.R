required <- c(
  "dplyr", "flextable", "ggplot2", "igraph", "knitr", "Matrix",
  "officedown", "officer", "reshape2", "rmarkdown", "RSpectra", "scales",
  "tidyr"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  install.packages(missing, dependencies = TRUE)
} else {
  message("All required R packages are already installed.")
}


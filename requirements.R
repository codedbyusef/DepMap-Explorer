# requirements.R
#
# Installs every package the DepMap Cell Death Explorer needs.
# Run with: Rscript requirements.R   (or source("requirements.R") from R)
#
# CRAN packages are installed only if missing. Bioconductor packages
# (depmap, ExperimentHub) go through BiocManager.

cran_pkgs <- c(
  "bslib",
  "shiny",
  "dplyr",
  "ggplot2",
  "DT",
  "shinycssloaders",
  "effsize"
)

bioc_pkgs <- c(
  "depmap",
  "ExperimentHub"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

still_missing <- c(cran_pkgs, bioc_pkgs)
still_missing <- still_missing[!vapply(still_missing, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}

message("All dependencies are installed.")

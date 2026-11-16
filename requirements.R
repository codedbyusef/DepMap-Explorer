# requirements.R
#
# Installs every package the DepMap Explorer needs.
# Run with: Rscript requirements.R   (or source("requirements.R") from R)
#
# CRAN packages are installed only if missing. Bioconductor packages go through
# BiocManager. GOfan is installed from GitHub.
#
# Runtime packages are needed to run the app. Build-only packages (depmap,
# ExperimentHub) are only used by preprocess.R, which writes the Parquet files
# in data/ that the app reads.

cran_pkgs <- c(
  "bslib",
  "shiny",
  "dplyr",
  "ggplot2",
  "DT",
  "shinycssloaders",
  "effsize",
  "arrow",
  "plotly",
  "remotes",
  "digest"  # checksums in get_data.R
)

bioc_pkgs <- c(
  "GO.db",
  "org.Hs.eg.db",
  # build-only (preprocess.R)
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

if (!requireNamespace("GOfan", quietly = TRUE)) {
  BiocManager::install("jianhong/GOfan", ask = FALSE, update = FALSE)
}

still_missing <- c(cran_pkgs, bioc_pkgs, "GOfan")
still_missing <- still_missing[!vapply(still_missing, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}

message("All dependencies are installed.")

# Download the preprocessed Parquet files from a GitHub Release into data/.
# Usage: Rscript get_data.R [release_tag]
# Files already present with a matching SHA-256 are skipped.
# To rebuild from DepMap instead, run `Rscript preprocess.R`.

REPO <- "codedbyusef/DepMap-Explorer"
TAG <- "data-24Q4"
DATA_DIR <- "data"
MANIFEST <- file.path(DATA_DIR, "checksums.txt")  # lines: "<sha256>  <file>"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) TAG <- args[1]

if (!file.exists(MANIFEST)) {
  stop("Missing ", MANIFEST, ". It is created when the release is made; see data/README.md.", call. = FALSE)
}
man <- read.table(MANIFEST, col.names = c("sha", "file"), colClasses = "character")

sha256 <- function(p) digest::digest(p, algo = "sha256", file = TRUE)
options(timeout = max(3600, getOption("timeout")))

for (i in seq_len(nrow(man))) {
  f <- man$file[i]
  dest <- file.path(DATA_DIR, f)
  if (file.exists(dest) && sha256(dest) == man$sha[i]) {
    message("ok       ", f)
    next
  }
  url <- sprintf("https://github.com/%s/releases/download/%s/%s", REPO, TAG, f)
  message("download ", f)
  tmp <- paste0(dest, ".part")
  utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
  if (sha256(tmp) != man$sha[i]) {
    unlink(tmp)
    stop("Checksum mismatch for ", f, call. = FALSE)
  }
  file.rename(tmp, dest)
}
message("All data files present and verified.")

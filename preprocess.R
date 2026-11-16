# =============================================================================
# preprocess.R
#
# One-off build step. Reads four DepMap Public release files (CSV, from Broad's
# portal or figshare), then writes compact Parquet files to data/ so the app
# never has to load them again. Run with:
#
#   Rscript preprocess.R <raw_dir>
#
# <raw_dir> holds CRISPRGeneEffect.csv, OmicsCNGene.csv,
# OmicsSomaticMutations.csv and Model.csv from one release (24Q4 or newer;
# the layout is the same in later releases). The Bioconductor depmap package
# is no longer used because it stops at 22Q2.
#
# To redo only the GO enrichment (for example after changing the cutoffs below),
# run:  Rscript preprocess.R go
# That step reads the existing lineage_effects.parquet and skips the slow
# DepMap tables.
#
# Outputs (all in data/, all gitignored):
#   crispr.parquet          depmap_id, gene_name, dependency   (sorted by gene)
#   meta.parquet            depmap_id, cell_line_name, primary_disease, lineage
#   copy_number.parquet     depmap_id, gene_name, log_copy_number (sorted by gene)
#   mutations.parquet       depmap_id, gene_name              (distinct, sorted)
#   lineage_effects.parquet gene_name, lineage, n_lineage, n_other, cohen_d, p, q
#   go_enrichment.parquet   lineage, ID, p.adjust, qvalue, count
#
# Sorting by gene_name and using small row groups lets arrow skip most of a
# file when the app filters to one gene.
#
# lineage_effects: Cohen's d is (lineage mean - other-lineage mean) / pooled SD,
# so a negative d means the lineage is MORE dependent than the rest. The p-value
# is a Welch t-test computed from summary statistics (fast for ~500k tests); it
# is not the Wilcoxon test used for the WT-vs-altered comparison in app.R. q is
# Benjamini-Hochberg across all genes within a lineage.
#
# go_enrichment: for each lineage, genes with d <= GO_D_CUT and q < GO_Q_CUT are
# tested for GO Biological Process enrichment (hypergeometric, universe = all
# screened genes, BH across all tested terms). Only terms in the cell death
# subtree (GO:0008219) are kept. qvalue = -log10(p.adjust), the scale GOfan
# expects for its fill.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
})

OUT_DIR <- "data"
MIN_LINEAGE_N <- 5
GO_D_CUT <- -0.5
GO_Q_CUT <- 0.05
GO_MIN_HITS <- 2
GO_MIN_TERM <- 5
GO_MAX_TERM <- 2000
CELL_DEATH_ROOT <- "GO:0008219"
CHUNK <- 100000

ARGS <- commandArgs(trailingOnly = TRUE)
GO_ONLY <- "go" %in% ARGS
RAW_DIR <- setdiff(ARGS, "go")[1]
if (!GO_ONLY && (is.na(RAW_DIR) || !dir.exists(RAW_DIR))) {
  stop("Usage: Rscript preprocess.R <raw_dir>   (or: Rscript preprocess.R go)", call. = FALSE)
}
raw <- function(f) {
  p <- file.path(RAW_DIR, f)
  if (!file.exists(p)) stop("Missing ", p, call. = FALSE)
  p
}

# Wide release matrix (models in rows, "SYMBOL (EntrezID)" in columns) to long
# depmap_id, gene_name, <value>. Drops NA values and duplicate symbols.
read_wide <- function(file, value) {
  dt <- data.table::fread(raw(file), check.names = FALSE)
  data.table::setnames(dt, 1, "depmap_id")
  long <- data.table::melt(dt, id.vars = "depmap_id", variable.name = "gene",
                           value.name = value, na.rm = TRUE)
  rm(dt)
  levels(long$gene) <- sub(" \\(.*$", "", levels(long$gene))
  long <- as_tibble(long) %>%
    mutate(gene_name = as.character(gene), .keep = "unused") %>%
    distinct(depmap_id, gene_name, .keep_all = TRUE) %>%
    relocate(depmap_id, gene_name)
  long
}

dir.create(OUT_DIR, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

write_sorted <- function(df, file) {
  df <- arrange(df, gene_name)
  write_parquet(df, out(file), compression = "zstd", chunk_size = CHUNK)
  message(sprintf("  %s: %s rows, %.1f MB", file,
                  format(nrow(df), big.mark = ","), file.size(out(file)) / 1e6))
}

if (!GO_ONLY) {
  # -- Raw tables ---------------------------------------------------------------
  message("CRISPR gene effects")
  crispr <- read_wide("CRISPRGeneEffect.csv", "dependency")
  write_sorted(crispr, "crispr.parquet")

  message("Cell line metadata")
  meta <- data.table::fread(raw("Model.csv"),
                            select = c("ModelID", "CellLineName", "OncotreePrimaryDisease", "OncotreeLineage")) %>%
    as_tibble() %>%
    transmute(depmap_id = ModelID, cell_line_name = CellLineName,
              primary_disease = OncotreePrimaryDisease, lineage = OncotreeLineage) %>%
    distinct(depmap_id, .keep_all = TRUE)
  write_parquet(meta, out("meta.parquet"), compression = "zstd")
  message(sprintf("  meta.parquet: %d rows", nrow(meta)))

  message("Mutation calls")
  # The release file lists coding variants only (moderate and high impact), so
  # any row counts as a mutation call.
  mut <- data.table::fread(raw("OmicsSomaticMutations.csv"), select = c("ModelID", "HugoSymbol")) %>%
    as_tibble() %>%
    transmute(depmap_id = ModelID, gene_name = HugoSymbol) %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    distinct()
  write_sorted(mut, "mutations.parquet")
  rm(mut)

  message("Copy number (large table)")
  # OmicsCNGene is log2(relative copy number + 1), diploid about 1.0, the same
  # scale as before, so CN_LOSS_CUT in app.R still applies.
  cn <- read_wide("OmicsCNGene.csv", "log_copy_number")
  write_sorted(cn, "copy_number.parquet")
  rm(cn)
  invisible(gc())

  # -- Lineage vs rest effect sizes -----------------------------------------------
  message("Lineage effect sizes")
  d <- crispr %>% inner_join(select(meta, depmap_id, lineage), by = "depmap_id") %>%
    filter(!is.na(lineage))
  rm(crispr)

  totals <- d %>%
    group_by(gene_name) %>%
    summarise(n_all = n(), s_all = sum(dependency), ss_all = sum(dependency^2), .groups = "drop")

  eff <- d %>%
    group_by(gene_name, lineage) %>%
    summarise(n_lineage = n(), s1 = sum(dependency), ss1 = sum(dependency^2), .groups = "drop") %>%
    filter(n_lineage >= MIN_LINEAGE_N) %>%
    inner_join(totals, by = "gene_name") %>%
    mutate(
      n_other = n_all - n_lineage,
      m1 = s1 / n_lineage,
      m2 = (s_all - s1) / n_other,
      v1 = pmax((ss1 - n_lineage * m1^2) / (n_lineage - 1), 0),
      v2 = pmax(((ss_all - ss1) - n_other * m2^2) / (n_other - 1), 0),
      pooled = sqrt(((n_lineage - 1) * v1 + (n_other - 1) * v2) / (n_lineage + n_other - 2)),
      cohen_d = (m1 - m2) / pooled,
      se = sqrt(v1 / n_lineage + v2 / n_other),
      tstat = (m1 - m2) / se,
      df = se^4 / ((v1 / n_lineage)^2 / (n_lineage - 1) + (v2 / n_other)^2 / (n_other - 1)),
      p = 2 * pt(-abs(tstat), df)
    ) %>%
    filter(is.finite(cohen_d), is.finite(p)) %>%
    group_by(lineage) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    select(gene_name, lineage, n_lineage, n_other, cohen_d, p, q) %>%
    arrange(lineage, gene_name)
  write_parquet(eff, out("lineage_effects.parquet"), compression = "zstd")
  message(sprintf("  lineage_effects.parquet: %s rows", format(nrow(eff), big.mark = ",")))

  # Spot check against the reference figure (kidney: VHL ~ 1.7, TP53 ~ 1.3).
  message("Kidney spot check (Cohen's d, q):")
  print(as.data.frame(
    eff %>% filter(lineage == "kidney", gene_name %in% c("VHL", "TP53", "MCL1", "PTEN", "GPX4", "BCL2L1")) %>%
      mutate(cohen_d = round(cohen_d, 2), q = signif(q, 2)) %>%
      select(gene_name, n_lineage, n_other, cohen_d, q)
  ))
} else {
  eff <- read_parquet(out("lineage_effects.parquet"))
}

# -- GO enrichment of selective dependencies ------------------------------------
message("GO enrichment")
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

universe <- unique(eff$gene_name)
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = universe, keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")
)) %>%
  filter(ONTOLOGYALL == "BP", !is.na(GOALL)) %>%
  distinct(SYMBOL, GOALL)

term_genes <- split(go_map$SYMBOL, go_map$GOALL)
term_size <- lengths(term_genes)
term_genes <- term_genes[term_size >= GO_MIN_TERM & term_size <= GO_MAX_TERM]
term_size <- lengths(term_genes)
death_terms <- c(CELL_DEATH_ROOT, GOBPOFFSPRING[[CELL_DEATH_ROOT]])
N <- length(universe)

enrich_one <- function(hits) {
  hits <- intersect(hits, universe)
  k <- length(hits)
  if (k < GO_MIN_HITS) return(NULL)
  x <- vapply(term_genes, function(g) sum(g %in% hits), numeric(1))
  keep <- x >= GO_MIN_HITS
  if (!any(keep)) return(NULL)
  p <- phyper(x[keep] - 1, term_size[keep], N - term_size[keep], k, lower.tail = FALSE)
  tibble(ID = names(x)[keep], count = as.integer(x[keep]), p = p)
}

go_all <- lapply(sort(unique(eff$lineage)), function(lin) {
  hits <- eff$gene_name[eff$lineage == lin & eff$cohen_d <= GO_D_CUT & eff$q < GO_Q_CUT]
  res <- enrich_one(hits)
  if (is.null(res)) return(NULL)
  res %>%
    mutate(lineage = lin, p.adjust = p.adjust(p, method = "BH")) %>%
    filter(ID %in% death_terms, p.adjust < 0.05) %>%
    mutate(qvalue = -log10(p.adjust)) %>%
    dplyr::select(lineage, ID, p.adjust, qvalue, count)
}) %>% bind_rows()

write_parquet(go_all, out("go_enrichment.parquet"), compression = "zstd")
message(sprintf("  go_enrichment.parquet: %d significant cell-death terms across %d lineages",
                nrow(go_all), n_distinct(go_all$lineage)))
message("Done.")

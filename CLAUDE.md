# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

- Install dependencies: `Rscript requirements.R` (CRAN packages, Bioconductor packages, and GOfan from GitHub).
- Download the prebuilt data files: `Rscript get_data.R` (fetches from the GitHub Release `data-22Q2` and checks `data/checksums.txt`).- Deploy to shinyapps.io: run `get_data.R` first so `data/*.parquet` are in the bundle (about 470 MB, under the 1 GB Free/Starter limit). The data are DepMap 22Q2, which is not the newest release; the Bioconductor `depmap` package stops there.
- Build the data files: `Rscript preprocess.R`. Run it once before the app (and again after a DepMap update). `Rscript preprocess.R go` redoes only the GO enrichment from the existing `lineage_effects.parquet`. It writes Parquet files to `data/` (gitignored). The app stops at startup with a message if they are missing.
- Run the app: `Rscript -e 'shiny::runApp(".", port = 4321, launch.browser = FALSE)'`.
- Syntax check after edits: `Rscript -e 'invisible(parse("app.R"))'`.
- There is no test suite and no linter configured.

The first `preprocess.R` run downloads the four DepMap tables into the ExperimentHub cache, so expect it to take a while. The app itself never calls depmap or ExperimentHub; it reads the Parquet files.

## Architecture

The app is `app.R` plus the one-off build script `preprocess.R`. `app.R` has four parts, in this order: constants, helper functions, the `ui` object (`page_sidebar` with a sidebar and a `navset_card_tab` of six tabs), and `server`.

**Data files.** `preprocess.R` writes `crispr`, `copy_number` and `mutations` Parquet files sorted by gene with small row groups, plus `meta`, `lineage_effects` and `go_enrichment`. `open_data()` in `app.R` runs once at startup into the global `DATA`, shared by all sessions. The three big tables stay on disk as arrow datasets; the small ones are read into memory. See `data/README.md`.

**Data flow.** `gene_frame()` takes the selected gene, filters the arrow datasets to that gene, `collect()`s only those rows, and returns one row per cell line. It joins CRISPR scores with metadata, copy number and mutation calls, then derives `altered`, `group` (WT/Altered) and `status`. `cohort()` filters `gene_frame()` by lineage. The first four tabs read from `cohort()` or `gene_frame()`; the last two read the precomputed tables. Keep this split: the gene filter is cached separately from the cancer filter, so changing the cancer type does not reload the gene.

**Altered definition.** A cell line is altered if it has a mutation call for the gene, or if `log_copy_number <= CN_LOSS_CUT`. The cutoff is 0.7311832. This is DepMap's own threshold for copy loss (heterozygous loss or deep deletion) on the log2(CN+1) scale, where diploid is about 1.0. Don't use a cutoff like 1.5 on this scale: it flags nearly every line as altered.

**Statistics.** `compare_groups()` runs `wilcox.test(alt, wt)` and `effsize::cohen.d(alt, wt)`. A positive d means WT lines are more dependent than altered lines. The sign was checked against group medians, so keep this direction.

**Dependency Plot tab.** Cohen's d for the genes picked in the tab, selected lineage vs all other lineages. d is lineage minus others, so negative means the lineage is more dependent. Significance is Wilcoxon with BH across the genes shown, computed from raw scores read from `crispr.parquet` for just those genes. `GENE_PANEL` (23 cell death genes in three roles taken from the project figure: death suppressor, executioner, growth suppressor / other) only defines the preset groups, all ticked by default. The figure's ferroptosis suppressor (AIFM2) and sensitizer (ACSL4) groups were dropped because each held one gene; add them back only with complete gene lists. Any other gene can be added and shows as "added gene". Ticking a group rewrites the selection as the ticked groups' genes plus added genes. The y-axis label colours use a vector `element_text()`, which ggplot2 warns is not officially supported. With the default groups, this reproduces the outlined genes in the project figure for kidney (VHL, TP53, MCL1, GPX4, BCL2L1, CASP9, PARP1) and the d values (VHL 1.69, TP53 1.28, GPX4 -0.67). The figure's test was inferred, not documented. `lineage_effects.parquet` is different: it uses a Welch test with BH across all genes and is only used to pick genes for GO enrichment.

**Cell Death GO tab.** `go_enrichment.parquet` holds, per lineage, GO Biological Process terms under cell death (GO:0008219) enriched among genes with d <= -0.5 and q < 0.05 (hypergeometric, BH). The tab draws them with GOfan (`getGraph`, `sunburstGO`; note the argument order is `df, org, g`), static ggplot2 or interactive plotly. The enrichment cutoffs are first guesses. Only 9 lineages have significant terms.

**Positive controls.** The project spec gives TP53 d = 1.83 and PTEN d = 0.63 as validation values. On pan-cancer data with the threshold above, the app gives TP53 d = -1.65 and PTEN d = -0.49. The spec values have not been reproduced, and their source is unknown. The header comment in `app.R` records this. Don't describe the controls as validated unless the numbers match.

**Defaults.** Cancer type is All Cancer Types, gene is BCL2L1, and the WT/altered overlay is on. These are set in `DEFAULT_CANCER` and `DEFAULT_GENE`, and the reset button restores them.

**UI behaviour that is deliberate.**
- The gene selectize uses client-side matching over the full gene list. `server = FALSE` and `maxOptions = 50` are set so partial matches work. Don't switch back to server-side mode, because it misses partial matches.
- The download buttons are inside `conditionalPanel(condition = "input.gene && input.gene !== ''")`, so they only appear when a gene is selected.
- `org.Hs.eg.db` is not attached, because it masks `dplyr::select`. Call it as `org.Hs.eg.db::org.Hs.eg.db`, and use `dplyr::select` after attaching it in scripts.
- The Dependency Plot and Cell Death GO tabs need a specific cancer type; with All Cancer Types they show a prompt.
- The Summary Statistics tab shows n, median, mean, min, max, the p-value, d and line counts. It has no interpretation sentence.

## Verification

Only the server logic and plot rendering (via `ggsave`) have been checked headlessly. The new tabs have not been checked in a browser. Source the app with `source("app.R", local = TRUE)`, then call `shiny::testServer(server, { session$setInputs(...) ... })`. `testServer` does not apply `update*` calls back to the inputs, so reset and default-setting behaviour needs a browser check.

## Repo

Default branch is `main`, remote `origin` is `github.com/codedbyusef/DepMap-Explorer` (public). `data/` holds only a README in git. The Parquet files are gitignored and built with `preprocess.R`.

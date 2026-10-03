# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

- Install dependencies: `Rscript requirements.R` (CRAN packages, plus depmap and ExperimentHub via BiocManager).
- Run the app: `Rscript -e 'shiny::runApp(".", port = 4321, launch.browser = FALSE)'`.
- Syntax check after edits: `Rscript -e 'invisible(parse("app.R"))'`.
- There is no test suite and no linter configured.

The first session after a fresh start downloads the four DepMap tables into the ExperimentHub cache. Expect a long first load. Later loads use the cache.

## Architecture

The whole app is `app.R`. It has four parts, in this order: constants, helper functions, the `ui` object (`page_sidebar` with a sidebar and a `navset_card_tab` of four tabs), and `server`.

**Data flow.** `datasets()` is a reactive that loads the four depmap tables once per session. The copy-number table is about 45M rows, so every session holds it in memory. `gene_frame()` takes the selected gene and returns one row per cell line. It joins CRISPR scores with metadata, copy number and mutation calls, then derives `altered`, `group` (WT/Altered) and `status`. `cohort()` filters `gene_frame()` by lineage. All four tabs read from `cohort()` or `gene_frame()`. Keep this split: the gene filter is cached separately from the cancer filter, so changing the cancer type does not reload the gene.

**Altered definition.** A cell line is altered if it has a mutation call for the gene, or if `log_copy_number <= CN_LOSS_CUT`. The cutoff is 0.7311832. This is DepMap's own threshold for copy loss (heterozygous loss or deep deletion) on the log2(CN+1) scale, where diploid is about 1.0. Don't use a cutoff like 1.5 on this scale: it flags nearly every line as altered.

**Statistics.** `compare_groups()` runs `wilcox.test(alt, wt)` and `effsize::cohen.d(alt, wt)`. A positive d means WT lines are more dependent than altered lines. The sign was checked against group medians, so keep this direction.

**Positive controls.** The project spec gives TP53 d = 1.83 and PTEN d = 0.63 as validation values. On pan-cancer data with the threshold above, the app gives TP53 d = -1.65 and PTEN d = -0.49. The spec values have not been reproduced, and their source is unknown. The header comment in `app.R` records this. Don't describe the controls as validated unless the numbers match.

**Defaults.** Cancer type is All Cancer Types, gene is BCL2L1, and the WT/altered overlay is on. These are set in `DEFAULT_CANCER` and `DEFAULT_GENE`, and the reset button restores them.

**UI behaviour that is deliberate.**
- The gene selectize uses client-side matching over the full gene list. `server = FALSE` and `maxOptions = 50` are set so partial matches work. Don't switch back to server-side mode, because it misses partial matches.
- The download buttons are inside `conditionalPanel(condition = "input.gene && input.gene !== ''")`, so they only appear when a gene is selected.
- The Summary Statistics tab shows n, median, mean, min, max, the p-value, d and line counts. It has no interpretation sentence.

## Verification

Only the server logic has been checked headlessly. Source the app with `source("app.R", local = TRUE)`, then call `shiny::testServer(server, { session$setInputs(...) ... })`. `testServer` does not apply `update*` calls back to the inputs, so reset and default-setting behaviour needs a browser check.

## Repo

Default branch is `main`, remote `origin` is `github.com/codedbyusef/shinyapps` (public). `data/` holds only a README. No data files are committed; everything loads through the depmap package.

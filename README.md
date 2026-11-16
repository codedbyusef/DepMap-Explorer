# DepMap Explorer

A Shiny app for exploring DepMap CRISPR gene dependency data. It compares cell lines with an altered gene (mutated or copy-number loss) against wild-type lines, and looks at how cell death regulators (death suppressors such as BCL2L1, GPX4 and MCL1, and executioners such as CASP3, BAX and MLKL) behave across cancer lineages.

The data are DepMap Public 24Q4.

## Features

Six tabs, driven by a cancer type and gene picked in the sidebar:

- Distribution: CRISPR scores for the selected cohort, with an optional WT/altered overlay and PNG download.
- Cell Line Table: per-line CRISPR score, copy number and mutation status, with CSV download.
- Summary Statistics: n, median, mean, min, max, Wilcoxon p-value and Cohen's d.
- Gene Summary: per-gene summary for the selected gene.
- Dependency Plot: Cohen's d for a panel of cell death genes, selected lineage against all other lineages, with BH-adjusted Wilcoxon significance.
- Cell Death GO: enriched GO Biological Process terms under cell death for each lineage, drawn with GOfan (static or interactive).

A cell line counts as altered if it has a mutation call for the gene, or if its log2 copy number is at or below DepMap's loss threshold (0.731). A positive Cohen's d means WT lines are more dependent than altered lines.

## Setup

```
Rscript requirements.R   # CRAN, Bioconductor and GOfan (from GitHub)
Rscript get_data.R       # downloads about 725 MB of Parquet files into data/
```

Instead of downloading, you can build the data yourself from the 24Q4 CSV files with `Rscript preprocess.R <raw_dir>`. It needs several GB of RAM, so run it under a memory cap (see [data/README.md](data/README.md)). See [data/README.md](data/README.md) for the file layout, source and license.

## Run

```
Rscript -e 'shiny::runApp(".", port = 4321, launch.browser = FALSE)'
```

Then open http://127.0.0.1:4321. The app stops at startup with a message if the Parquet files are missing.

## How data is stored

The three large tables (`crispr`, `copy_number`, `mutations`) are Parquet files sorted by gene with small row groups. The app opens them as arrow datasets and reads only the rows for the selected gene. The small tables (`meta`, `lineage_effects`, `go_enrichment`) are loaded into memory at startup. The app never calls `depmap` or ExperimentHub at runtime.

## Deploying to shinyapps.io

Run `Rscript get_data.R` first so the Parquet files go into the bundle, then `rsconnect::deployApp()`. The bundle is about 725 MB, under the 1 GB Free/Starter limit.

## Known limitations

- The project spec gives TP53 d = 1.83 and PTEN d = 0.63 as positive controls. The app gives TP53 d = -1.76 and PTEN d = -0.40 on 24Q4 (-1.65 and -0.49 on 22Q2). I tried mutation-only, mutation or copy loss at several cutoffs, and copy loss only, on 22Q2 and on 24Q4. TP53 stays between -1.55 and -1.76 and PTEN between -0.06 and -0.71. The sign differs from the spec (the app uses altered minus WT) and no definition reaches 1.83, so the controls are not validated and the spec's source is unknown.
- The GO enrichment cutoffs are an analysis choice with no single right value. With q < 0.05, d <= -0.3, -0.5 and -0.8 gave significant cell death terms in 10, 9 and 13 of 25 lineages on 22Q2. Stricter cutoffs find more lineages because the gene sets get smaller and more specific. The shipped data use d <= -0.5.

## Data license

The data are derived from DepMap (Broad Institute, https://depmap.org) under CC BY 4.0, with changes described in [data/README.md](data/README.md).

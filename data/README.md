# Data

The DepMap Cell Death Explorer does not ship any data files. Everything is loaded at runtime from **Bioconductor** through **ExperimentHub**, using the `depmap` package. No manual downloads and no API tokens are needed.

## Datasets

On startup, `app.R` calls these four functions once per session:

| Function | Contents |
|---|---|
| `depmap::depmap_crispr()` | CRISPR gene effect scores (Chronos) per cell line and gene |
| `depmap::depmap_copyNumber()` | Copy number per cell line and gene |
| `depmap::depmap_mutationCalls()` | Mutation calls per cell line and gene |
| `depmap::depmap_metadata()` | Cell line metadata, including lineage and primary disease |

## First launch

The first run downloads the datasets into the ExperimentHub cache. They are large, so this can take a while and requires an internet connection. A spinner shows while the data loads. Later runs read from the local cache.

## Requirements

Run `Rscript requirements.R` from the project root to install the packages the app needs.

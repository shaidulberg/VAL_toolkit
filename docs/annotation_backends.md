# Optional annotation backends

The ACS workflow is designed around multiple annotation methods. In the manuscript-style workflow, each method contributes one vote for a broad immune lineage, and ACS ranking is then computed from:

1. VAL: the number of methods assigning the consensus/majority lineage.
2. Mean consensus confidence: the mean confidence among methods that assigned that consensus lineage.

This repository supports two ways to provide annotation methods.

## Recommended stable mode: annotation columns

Users can provide annotation labels and confidence scores that already exist in `adata.obs`, or as an external CSV. This is the most stable way to compute manuscript-style ACS today.

Example columns:

```text
celltypist_label
celltypist_confidence
singler_hpca_label
singler_hpca_confidence
singler_monaco_label
singler_monaco_confidence
scimilarity_label
scimilarity_confidence
```

These columns are then listed in `configs/acs_figures.example.yaml` under `annotation_methods`.

## Backend runner

The backend runner is:

```bash
python scripts/run_annotation_backends.py --config configs/annotation_backends.example.yaml
```

It writes:

```text
annotation_backend_columns.csv
annotation_backend_status.csv
annotation_backends_annotated.h5ad
```

The completed annotation columns can then be used by:

```bash
python scripts/run_acs_figures.py --config configs/acs_figures.example.yaml
```

## Current backend status

### CellTypist

Implemented as an optional Python backend.

The backend reads an `.h5ad`, runs CellTypist, and writes:

```text
celltypist_label
celltypist_confidence
```

`celltypist_confidence` is the CellTypist probability assigned to the selected label when available. If the selected label is not found in the probability matrix, the row maximum probability is used as a fallback.

### SingleR

Implemented as an optional R/Bioconductor backend.

The backend exports the selected AnnData expression matrix to MatrixMarket format, runs SingleR through `Rscript`, and writes one label/confidence pair per configured celldex reference. See `docs/singler_backend.md` for R package installation and configuration details.

### SCimilarity

Implemented as an optional Python backend.

The backend uses the SCimilarity `CellAnnotation` API with a user-provided local model directory. It prepares the selected AnnData matrix, aligns the query genes to the model gene order, computes SCimilarity embeddings, obtains KNN-based cell-type predictions, and writes:

```text
scimilarity_label
scimilarity_confidence
```

`scimilarity_confidence` defaults to the SCimilarity `vsAll` statistic. Users can also request `vs2nd`, `vsAll_weighted`, or `vs2nd_weighted` through the config.

SCimilarity remains optional because it requires a separate package install and downloaded model files. If either is missing, set `on_unavailable: skip` to continue with the other completed methods. See `docs/scimilarity_backend.md`.

## Why all backends are optional

The core ACS workflow should remain usable even when optional annotation packages are unavailable. The runner therefore records a status for each backend and writes the methods that completed successfully. This prevents a missing R or SCimilarity installation from blocking CellTypist or ACS workflows.

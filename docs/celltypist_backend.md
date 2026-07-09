# Optional CellTypist backend

This backend is the first automatic annotation backend for the VAL/ACS toolkit.
It is optional: the core ACS workflow can still run from annotation columns that
are already present in `adata.obs` or in an external annotation CSV.

## What this backend does

The CellTypist backend reads one `.h5ad` file, runs CellTypist, and writes:

- `celltypist_label`
- `celltypist_confidence`

These columns can then be used as one annotation method in the manuscript-style
ACS workflow.

## Basic command

```bash
python scripts/create_example_h5ad.py
python scripts/run_celltypist_annotations.py --config configs/celltypist_annotations.example.yaml
```

This produces:

```text
figures/example_celltypist_annotations/celltypist_annotations.csv
figures/example_celltypist_annotations/celltypist_annotated.h5ad
```

## Input expectations

CellTypist expects a cell-by-gene matrix with gene symbols in `adata.var_names`.
For standard built-in models, expression should be log1p-normalized to 10,000
counts per cell. If your `.h5ad` already contains that matrix in `adata.X`, leave
preprocessing disabled. If your `.h5ad` contains raw counts, use:

```yaml
expression:
  normalize_to_target_sum: true
  target_sum: 10000
  log1p: true
```

## How this connects to ACS

After running CellTypist, use the resulting `celltypist_annotated.h5ad` as input
to an ACS config, together with any other annotation methods available in
`adata.obs` or an annotation CSV.

A single CellTypist annotation alone is not sufficient to reproduce the full
manuscript VAL/ACS logic because VAL is defined as cross-method agreement. It is
one vote in the multi-method ACS framework.

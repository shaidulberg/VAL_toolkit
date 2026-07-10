# Automatic SingleR backend

The SingleR backend is optional because it requires an R/Bioconductor installation in addition to the Python environment.

## What it does

When enabled, the backend:

1. reads the input `.h5ad`,
2. selects the configured expression matrix,
3. optionally normalizes each cell to 10,000 counts and applies `log1p`,
4. exports the query matrix to MatrixMarket format,
5. runs SingleR through `Rscript` for each configured celldex reference,
6. writes one label/confidence pair per reference into the backend-annotated `.h5ad`, and
7. passes those columns into the ACS workflow.

Example output columns:

```text
singler_hpca_label
singler_hpca_confidence
singler_monaco_label
singler_monaco_confidence
singler_dice_label
singler_dice_confidence
singler_encode_label
singler_encode_confidence
singler_hema_label
singler_hema_confidence
```

The default confidence value is the SingleR score for the assigned label. Set `score_column: delta.next` to use the margin to the next label instead.

## R dependencies

Install the required Bioconductor packages in R:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "SingleR",
  "celldex",
  "Matrix",
  "SummarizedExperiment"
))
```

Then confirm that the same `Rscript` visible to the terminal can load the packages:

```bash
Rscript -e "cat(all(sapply(c('SingleR','celldex','Matrix','SummarizedExperiment'), requireNamespace, quietly=TRUE)))"
```

The command should print `TRUE`.

## Enable SingleR

In `configs/acs_pipeline.example.yaml` or `configs/annotation_backends.example.yaml`:

```yaml
annotation_backends:
  singler:
    enabled: true
    references:
      - HPCA
      - Monaco
      - DICE
      - BlueprintEncode
      - Novershtern
    prefix: singler
    rscript: Rscript
    label_field: label.main
    prune_labels: true
    score_column: assigned_label_score
    assay_type_ref: logcounts
    keep_temp: false
    on_unavailable: skip
```

Supported reference aliases:

| Config name | celldex reference |
|---|---|
| `HPCA` | `HumanPrimaryCellAtlasData()` |
| `Monaco` | `MonacoImmuneData()` |
| `DICE` | `DatabaseImmuneCellExpressionData()` |
| `BlueprintEncode` or `ENCODE` | `BlueprintEncodeData()` |
| `Novershtern` or `HEMA` | `NovershternHematopoieticData()` |

## Run

```bash
python scripts/run_acs_pipeline.py --config configs/acs_pipeline.example.yaml
```

If `pipeline.auto_include_backend_methods: true`, completed SingleR columns are automatically added to the generated ACS config. Users do not need to manually list every `singler_*_label` column.

## Notes

- SingleR can be slow on large datasets. Start with one or two references while testing.
- The current bridge exports a MatrixMarket file as an interoperability layer. This is robust but can create large temporary files.
- Use `keep_temp: true` only for debugging; otherwise temporary matrix files are removed after the run.

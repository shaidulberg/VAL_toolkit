# ACS config schema

This page summarizes the main fields used by `configs/acs_figures.example.yaml`.

## Required fields

```yaml
input_h5ad: /path/to/data.h5ad
output_dir: figures/my_acs_run
```

## Metadata columns

```yaml
obs_columns:
  patient_id: patient_id
  response: response
  timepoint: timepoint
  sample_id:
```

`patient_id` and `response` are required. `timepoint` and `sample_id` are optional. If `sample_id` is omitted and `timepoint` is present, samples are defined as `patient_id + timepoint`.

## Response groups

```yaml
response_groups:
  positive: R
  negative: NR
```

Values must match the contents of the response column.

## Major cell types

```yaml
major_cell_types:
  - B
  - CD4 T
  - CD8 T
  - DC
  - Monocytes
  - NK
```

These are the supported ACS figure panels.

## Annotation source

### Mode 1: annotation columns in adata.obs

```yaml
annotation_source:
  mode: obs_columns
```

### Mode 2: external annotation CSV

```yaml
annotation_source:
  mode: annotation_csv
  path: /path/to/annotation_calls.csv
  cell_id_column: cell_id
```

### Mode 3: self-contained marker fallback

```yaml
annotation_source:
  mode: marker_ensemble
```

This mode does not require external annotation methods, but it is ACS-like rather than the exact multi-method manuscript ACS calculation.

## Annotation methods

```yaml
annotation_methods:
  - name: CellTypist
    label_column: celltypist_label
    confidence_column: celltypist_confidence

  - name: SingleR_HPCA
    label_column: singler_hpca_label
    confidence_column: singler_hpca_confidence
```

`confidence_column` is optional. If omitted, confidence is treated as 1.0.

## Label mapping

Grouped style:

```yaml
label_mapping:
  B:
    - B cell
    - B cells
    - Memory B
  CD4 T:
    - CD4 T cell
    - Treg
```

Direct style:

```yaml
label_mapping:
  B cells: B
  CD4 T cells: CD4 T
```

The grouped style is recommended for readability.

## Binning

```yaml
binning:
  n_bins: 10
  denominator: all_cells
```

`denominator` can be `all_cells` or `cell_type_cells`.

## Statistics

```yaml
statistics:
  fdr_scope: global
```

`fdr_scope` can be `global` or `within_cell_type`.

## Plot settings

```yaml
plot:
  fdr_threshold: 0.05
  dpi: 300
  basename: acs_figures_A_to_F
```

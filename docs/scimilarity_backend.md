# SCimilarity backend

The SCimilarity backend is optional and uses the Python SCimilarity API. It is intended for users who already have:

1. the `scimilarity` Python package installed, and
2. a downloaded local SCimilarity model directory.

The backend is disabled by default because the model files are external and can be large.

## Configuration

In `configs/acs_pipeline.example.yaml` or `configs/annotation_backends.example.yaml`:

```yaml
annotation_backends:
  scimilarity:
    enabled: true
    model_path: /path/to/scimilarity/model
    prefix: scimilarity
    use_gpu: false
    k: 50
    ef: 100
    weighting: false
    confidence_column: vsAll
    gene_overlap_threshold: 5000
    safelist: []
    blocklist: []
    on_unavailable: skip
```

## Output columns

When completed, the backend writes:

```text
scimilarity_label
scimilarity_confidence
```

It also writes optional diagnostic columns when available, such as:

```text
scimilarity_min_dist
scimilarity_max_dist
scimilarity_vs2nd
scimilarity_vsAll
scimilarity_vs2nd_weighted
scimilarity_vsAll_weighted
```

The one-command ACS pipeline automatically includes completed backend columns when `pipeline.auto_include_backend_methods: true`.

## Expression preprocessing

SCimilarity expects log-normalized expression aligned to the model gene order. The backend uses the shared `expression:` section to select and optionally normalize/log-transform the matrix before alignment. For raw count-like matrices, use:

```yaml
expression:
  normalize_to_target_sum: true
  target_sum: 10000
  log1p: true
```

For an already log-normalized matrix, use:

```yaml
expression:
  normalize_to_target_sum: false
  log1p: false
```

## Optional safelist/blocklist

The backend exposes SCimilarity safelist and blocklist options. These are useful when restricting predictions to immune labels or excluding broad labels. Example:

```yaml
safelist:
  - B cell
  - CD4-positive, alpha-beta T cell
  - CD8-positive, alpha-beta T cell
```

Leave both empty for the default full SCimilarity label space.

## Failure behavior

Set:

```yaml
on_unavailable: skip
```

to continue with CellTypist/SingleR if SCimilarity is unavailable. Set `on_unavailable: error` when SCimilarity is required.

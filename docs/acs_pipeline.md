# One-command ACS pipeline

The one-command pipeline runs the two main steps needed for ACS figures:

1. Optional annotation backends, currently automatic CellTypist and scaffolded SingleR/SCimilarity.
2. ACS-ranked bin analysis and ACS figures A-F.

```bash
python scripts/run_acs_pipeline.py --config configs/acs_pipeline.example.yaml
```

The command writes three output folders under `output_dir`:

- `annotation_backends/`: backend-specific annotation columns, status table, and optionally an annotated `.h5ad`.
- `acs_figures/`: ACS cell table, ACS-ranked bins, abundance table, association results, and figures A-F.
- `generated_configs/`: the backend and ACS configs generated from the single pipeline config.

## Required user inputs

At minimum, the user must provide:

- an `.h5ad` file,
- an `obs` column identifying patient or sample,
- an `obs` response column if responder/non-responder ACS association figures are desired.

The default example enables CellTypist and then runs ACS using CellTypist's label/confidence columns.
This is a useful one-file automatic workflow, but with only one annotation method, `VAL` can only be 1.
For manuscript-style multi-method ACS, provide additional annotation columns from SingleR, SCimilarity, or other methods and list them under `annotation_methods` in the config.

## Skipping annotation backends

If the `.h5ad` already contains annotation label/confidence columns, run:

```bash
python scripts/run_acs_pipeline.py --config configs/acs_pipeline.example.yaml --skip-backends
```

In this mode, ACS uses the existing annotation columns directly.

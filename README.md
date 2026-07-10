# VAL / ACS single-cell toolkit

This repository contains a user-facing implementation of **ACS-ranked single-cell response-association figures** for immune single-cell datasets stored as `.h5ad` files.

In the manuscript, cells are first evaluated by **VAL**: the number of annotation methods assigning the same majority/consensus label. Cells are then ranked within each VAL stratum by the mean confidence score among the methods that support that consensus label. This two-level weighted ranking is referred to here as the **Aggregated Confidence Score (ACS)** ranking.

The main public workflow is:

```text
user h5ad
  + patient/sample metadata
  + response labels
  + annotation-method labels/confidences in adata.obs or an external CSV
  → harmonize each method's labels to B/CD4 T/CD8 T/DC/Monocytes/NK/Other
  → compute consensus cell type and VAL
  → rank cells by ACS within each major lineage
  → split each lineage into non-overlapping 10% ACS-ranked bins
  → compare ranked-bin abundance between response groups
  → generate ACS figures A-F for the user's own dataset
```

The repository also includes a self-contained marker-ensemble fallback for users who do not yet have CellTypist/SingleR/SCimilarity-style annotation columns. That fallback can generate ACS-like figures, but the manuscript-faithful ACS workflow uses multiple annotation-method columns.

---

## Repository structure

```text
.
├── configs/
│   ├── acs_figures.example.yaml
│   ├── val_ranked_bin_panels.example.yaml
│   ├── annotation_benchmark.example.yaml
│   └── auc_barplot.example.yaml
├── docs/
│   ├── acs_algorithm.md
│   └── acs_config_schema.md
├── example_data/
│   ├── annotation_predictions_example.csv
│   └── auc_barplot_input_example.csv
├── scripts/
│   ├── create_example_h5ad.py
│   ├── run_acs_figures.py
│   ├── run_val_ranked_bin_panels.py
│   ├── run_annotation_benchmark.py
│   └── run_auc_barplot.py
├── src/val_toolkit/
│   ├── acs.py
│   ├── annotation_harmonization.py
│   ├── val_panels.py
│   ├── annotation_metrics.py
│   ├── config.py
│   ├── labels.py
│   ├── plotting.py
│   └── validation.py
├── tests/
├── environment.yml
├── pyproject.toml
└── requirements.txt
```

---

## Installation

### Option 1: conda

```bash
conda env create -f environment.yml
conda activate val-single-cell
pip install -e .
```

### Option 2: pip / venv

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

---


## Recommended one-command workflow

For most users, start with the one-command ACS pipeline. It runs enabled annotation backends and then generates ACS figures A-F.

```bash
conda env create -f environment.yml
conda activate val-single-cell

python scripts/create_example_h5ad.py
python scripts/run_acs_pipeline.py --config configs/acs_pipeline.example.yaml
```

For a user's own dataset, edit `configs/acs_pipeline.example.yaml` to point `input_h5ad` to the user's `.h5ad` file and set the metadata columns under `obs_columns`.

The current automatic backends are CellTypist and SingleR. SingleR runs through an optional R/Bioconductor bridge; SCimilarity is scaffolded and can be represented through precomputed annotation columns until its automatic runner is implemented. With only CellTypist enabled, the pipeline ranks cells by CellTypist confidence; manuscript-style multi-method ACS requires additional annotation-method columns.

See `docs/acs_pipeline.md` for details.

## Quick test: generate ACS figures A-F from an example h5ad

Create a tiny synthetic `.h5ad` file. The example includes expression data, sample metadata, response labels, and fake multi-method annotation columns so the ACS workflow can be tested without installing external annotation packages.

```bash
python scripts/create_example_h5ad.py
```

Run the ACS workflow:

```bash
python scripts/run_acs_figures.py --config configs/acs_figures.example.yaml
```

Outputs are written to:

```text
figures/example_acs_figures/
```

Key outputs:

```text
acs_cell_table.csv
acs_ranked_cell_bins.csv
acs_ranked_bin_abundance_by_sample.csv
acs_ranked_bin_association_results.csv
acs_figures_A_to_F.png/pdf/svg
ranked_bin_panel_A_B.png/pdf/svg
ranked_bin_panel_B_CD4_T.png/pdf/svg
...
```

Run tests:

```bash
pytest -q
```

---

## Run ACS figures A-F on your own h5ad

Your `.h5ad` must contain sample-level metadata in `adata.obs`.

At minimum, the ACS figure workflow needs:

| Required item | Example column | Why it is needed |
|---|---:|---|
| patient/sample ID | `patient_id` | aggregates cells by sample |
| response group | `response` | compares responders versus non-responders |
| optional timepoint | `timepoint` | separates pre/post or other repeated samples |

For manuscript-faithful ACS, your `.h5ad` should also contain annotation labels and confidence scores from multiple methods in `adata.obs`, for example:

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

Copy and edit the config:

```bash
cp configs/acs_figures.example.yaml configs/acs_figures.mydata.yaml
```

Edit these fields:

```yaml
input_h5ad: /path/to/my_dataset.h5ad
output_dir: figures/my_dataset_acs_figures

obs_columns:
  patient_id: patient_id
  timepoint: timepoint
  response: response

response_groups:
  positive: R
  negative: NR
```

Edit `annotation_methods` to match your annotation columns:

```yaml
annotation_methods:
  - name: CellTypist
    label_column: celltypist_label
    confidence_column: celltypist_confidence

  - name: SingleR_HPCA
    label_column: singler_hpca_label
    confidence_column: singler_hpca_confidence

  - name: SCimilarity
    label_column: scimilarity_label
    confidence_column: scimilarity_confidence
```

Then run:

```bash
python scripts/run_acs_figures.py --config configs/acs_figures.mydata.yaml
```

You can also override common settings from the command line:

```bash
python scripts/run_acs_figures.py \
  --config configs/acs_figures.mydata.yaml \
  --h5ad /path/to/my_dataset.h5ad \
  --patient-col patient_id \
  --timepoint-col timepoint \
  --response-col response \
  --output-dir figures/my_dataset_acs_figures
```

---

## ACS definition used here

For each cell:

1. Harmonize each annotation method's raw label to a major immune lineage: `B`, `CD4 T`, `CD8 T`, `DC`, `Monocytes`, `NK`, or `Other`.
2. Determine the consensus label as the majority harmonized label, excluding `Other` unless all methods are `Other`.
3. Define `VAL` as the number of annotation methods that assign the consensus label.
4. Compute `mean_consensus_confidence` as the mean confidence score across only the methods that assign the consensus label.
5. Rank cells within each major cell type by `VAL` descending, then by `mean_consensus_confidence` descending.
6. Divide cells in each major cell type into non-overlapping 10% ACS-ranked bins.
7. For each bin, compute per-sample abundance and compare response groups using two-sided Mann-Whitney U tests.
8. Report FDR-adjusted q-values and plot `-log10(FDR q-value)` across ACS-ranked bins.

More detail is provided in [`docs/acs_algorithm.md`](docs/acs_algorithm.md).

---

## Annotation CSV mode

If your annotation calls are not stored in `adata.obs`, you can provide an external CSV matched by cell ID:

```yaml
annotation_source:
  mode: annotation_csv
  path: /path/to/annotation_calls.csv
  cell_id_column: cell_id
```

The CSV must contain one row per cell and the label/confidence columns listed under `annotation_methods`.

---

## Self-contained marker-ensemble fallback

If you only have expression data and metadata but no annotation-method outputs, use:

```yaml
annotation_source:
  mode: marker_ensemble
```

This mode computes marker-based votes for broad immune labels and generates ACS-like ranked-bin figures. It is useful for quick testing and fully self-contained runs, but it is not the exact multi-method manuscript ACS calculation.

---

## Additional helper workflows

### Annotation benchmark

Prepare a CSV where each row is a cell and columns include a ground-truth label plus one prediction column per method. Then run:

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.example.yaml
```

### AUC comparison bar plot

Prepare a CSV with at least these columns:

```text
Dataset, Signature, AUC
```

Then run:

```bash
python scripts/run_auc_barplot.py --config configs/auc_barplot.example.yaml
```

The AUC bar plot is a manuscript-summary workflow. It is not the default bring-your-own-h5ad workflow because reproducing the full external-validation bar plot requires multiple validation datasets and response labels.

## Optional annotation backends

The toolkit includes an optional backend runner that is designed for the full manuscript-style annotation framework. CellTypist is implemented as an automatic Python backend. SingleR is implemented through an optional R/Bioconductor bridge. SCimilarity is scaffolded and documented, but is not yet automatically executed because it requires additional model/environment setup.

```bash
python scripts/create_example_h5ad.py
python scripts/run_annotation_backends.py --config configs/annotation_backends.example.yaml
```

The backend runner writes completed annotation columns, a backend status table, and optionally an annotated h5ad. Completed columns can then be used by the ACS workflow.

```bash
python scripts/run_acs_figures.py --config configs/acs_figures.example.yaml
```

Current backend status:

- **CellTypist**: implemented.
- **SingleR**: implemented through optional `Rscript` + Bioconductor packages; see `docs/singler_backend.md`.
- **SCimilarity**: scaffolded; provide precomputed SCimilarity labels/confidences through `adata.obs` or an annotation CSV for now.

See `docs/annotation_backends.md and docs/singler_backend.md` for details.

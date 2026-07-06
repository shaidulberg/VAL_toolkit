# VAL single-cell annotation and response-signature analysis

This repository contains runnable code for reproducing core analyses from the VAL single-cell annotation workflow, including:

- one-vs-rest annotation benchmark metrics across cell types and methods;
- publication-style benchmark bar plots with mean ± SEM and per-dataset points;
- response-signature AUC comparison bar plots;
- configurable input paths through YAML files, with no hard-coded local paths.

The repository is designed so another user can run the analyses by editing config files rather than editing Python source code.

---

## Repository structure

```text
.
├── configs/
│   ├── annotation_benchmark.example.yaml
│   └── auc_barplot.example.yaml
├── example_data/
│   ├── annotation_predictions_example.csv
│   └── auc_barplot_input_example.csv
├── scripts/
│   ├── run_annotation_benchmark.py
│   └── run_auc_barplot.py
├── src/
│   └── val_toolkit/
│       ├── annotation_metrics.py
│       ├── config.py
│       ├── labels.py
│       ├── plotting.py
│       └── validation.py
├── tests/
│   └── test_annotation_metrics.py
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

## Quick test with example data

Run the annotation benchmark example:

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.example.yaml
```

Run the AUC bar plot example:

```bash
python scripts/run_auc_barplot.py --config configs/auc_barplot.example.yaml
```

Example outputs are written to:

```text
figures/example_annotation_benchmark/
figures/example_auc_barplot/
```

---

## Running on your own data

### Annotation benchmark

Prepare a CSV where each row is a cell and columns include:

- a ground-truth label column, for example `Reference`;
- one column per annotation method, for example `CellTypist`, `SCimilarity`, `singler.ENCODE`.

Then copy and edit:

```bash
cp configs/annotation_benchmark.example.yaml configs/annotation_benchmark.mydata.yaml
```

Update these fields:

```yaml
datasets:
  - name: HN
    annotation_csv: /path/to/your/annotation_table.csv
    truth_col: Reference
    method_cols:
      - CellTypist
      - SCimilarity
      - singler.ENCODE

cell_types:
  - B
  - CD4
  - CD8
  - DC
  - Mono
  - NK
```

Then run:

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.mydata.yaml
```

### AUC comparison bar plot

Prepare a CSV with at least these columns:

```text
Dataset, Signature, AUC
```

Then copy and edit:

```bash
cp configs/auc_barplot.example.yaml configs/auc_barplot.mydata.yaml
```

Run:

```bash
python scripts/run_auc_barplot.py --config configs/auc_barplot.mydata.yaml
```

---

## Input notes

- The plotting code uses raw AUC values as provided in the input CSV.
- Error bars are SEM across datasets.
- Individual dataset points are overlaid on each bar.
- Display labels are cleaned automatically, for example `singler.ENCODE` becomes `SingleR.ENCODE`, `CD4` becomes `CD4 T`, and `Mono` becomes `Monocytes`.
- The code does not discover files automatically. All input paths are explicitly specified in YAML config files.

---

## Reproducibility notes

For a manuscript release, recommended additions before public upload are:

1. commit the exact config files used for the final figures;
2. include a `data_availability.md` file describing where raw public datasets can be downloaded;
3. include non-sensitive derived input CSVs if redistribution is permitted;
4. add the final figure-generating commands to this README.

---

## Citation

If you use this code, please cite the associated manuscript once available.

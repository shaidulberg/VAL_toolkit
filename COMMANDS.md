# Common commands

## Create and activate environment

```bash
conda env create -f environment.yml
conda activate val-single-cell
pip install -e .
```

## Test the ACS h5ad workflow

```bash
python scripts/create_example_h5ad.py
python scripts/run_acs_figures.py --config configs/acs_figures.example.yaml
pytest -q
```

## Run ACS figures on your own h5ad

```bash
cp configs/acs_figures.example.yaml configs/acs_figures.mydata.yaml
# edit configs/acs_figures.mydata.yaml
python scripts/run_acs_figures.py --config configs/acs_figures.mydata.yaml
```

## Marker-ensemble fallback

Set this in `configs/acs_figures.mydata.yaml`:

```yaml
annotation_source:
  mode: marker_ensemble
```

Then run:

```bash
python scripts/run_acs_figures.py --config configs/acs_figures.mydata.yaml
```

## Legacy VAL-ranked bin workflow

```bash
python scripts/create_example_h5ad.py
python scripts/run_val_ranked_bin_panels.py --config configs/val_ranked_bin_panels.example.yaml
```

## Annotation benchmark helper

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.example.yaml
```

## AUC barplot helper

```bash
python scripts/run_auc_barplot.py --config configs/auc_barplot.example.yaml
```

## Git update workflow

```bash
git status
git add .
git commit -m "Describe what changed"
git push
```

## Optional annotation backends

```bash
conda activate val-single-cell
python scripts/create_example_h5ad.py
python scripts/run_annotation_backends.py --config configs/annotation_backends.example.yaml
```

CellTypist is implemented. SingleR and SCimilarity are scaffolded and should be supplied as precomputed annotation columns or CSVs until their automatic backends are implemented.

## Legacy CellTypist-only runner

```bash
python scripts/run_celltypist_annotations.py --config configs/celltypist_annotations.example.yaml
```

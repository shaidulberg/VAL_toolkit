# Reproducibility commands

Run these from the repository root after installing the environment.

## Example annotation benchmark

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.example.yaml
```

## Example AUC bar plot

```bash
python scripts/run_auc_barplot.py --config configs/auc_barplot.example.yaml
```

## Tests

```bash
pytest -q
```

## Manuscript runs

Add the exact final manuscript commands here, for example:

```bash
python scripts/run_annotation_benchmark.py --config configs/annotation_benchmark.manuscript.yaml
python scripts/run_auc_barplot.py --config configs/B_40_50_auc.manuscript.yaml
```

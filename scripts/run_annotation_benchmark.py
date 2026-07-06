#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pandas as pd

from val_toolkit.annotation_metrics import DatasetConfig, compute_metrics_from_dataset_configs
from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.plotting import plot_annotation_metric_by_cell_type, plot_annotation_metric_by_method
from val_toolkit.validation import ensure_output_dir, read_csv_checked


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute and plot one-vs-rest annotation benchmark metrics."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    output_dir = ensure_output_dir(resolve_path(config.get("output_dir", "figures/annotation_benchmark"), base_dir))
    metrics_to_plot = config.get("metrics_to_plot", ["F1", "Precision", "Recall"])
    ylim = tuple(config.get("ylim", [0, 1.05]))
    dpi = int(config.get("dpi", 300))
    cell_types = config.get("cell_types")
    if not cell_types:
        raise ValueError("Config must define cell_types.")

    method_order = config.get("method_order")
    cell_type_order = config.get("cell_type_order", cell_types)

    if "metrics_csv" in config and config["metrics_csv"]:
        metrics_path = resolve_path(config["metrics_csv"], base_dir)
        metrics = read_csv_checked(metrics_path, required_columns=["Dataset", "Method", "Cell Type", "Metric", "Value"])
    else:
        dataset_configs: list[DatasetConfig] = []
        for dataset in config.get("datasets", []):
            dataset_configs.append(
                DatasetConfig(
                    name=dataset["name"],
                    annotation_csv=resolve_path(dataset["annotation_csv"], base_dir),
                    truth_col=dataset["truth_col"],
                    method_cols=list(dataset["method_cols"]),
                )
            )
        if not dataset_configs:
            raise ValueError("Config must define either metrics_csv or datasets.")

        metrics = compute_metrics_from_dataset_configs(
            datasets=dataset_configs,
            cell_types=cell_types,
            drop_missing_predictions=bool(config.get("drop_missing_predictions", True)),
        )
        metrics.to_csv(output_dir / "annotation_metrics_raw_precision_recall_f1_accuracy.csv", index=False)

    for metric in metrics_to_plot:
        plot_annotation_metric_by_method(
            metrics=metrics,
            output_dir=output_dir,
            metric=metric,
            cell_type_order=cell_type_order,
            method_order=method_order,
            ylim=ylim,
            dpi=dpi,
        )
        plot_annotation_metric_by_cell_type(
            metrics=metrics,
            output_dir=output_dir,
            metric=metric,
            cell_type_order=cell_type_order,
            method_order=method_order,
            ylim=ylim,
            dpi=dpi,
        )

    print(f"Saved annotation benchmark outputs to: {output_dir}")


if __name__ == "__main__":
    main()

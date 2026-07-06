#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.plotting import plot_auc_barplot
from val_toolkit.validation import ensure_output_dir, read_csv_checked


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot response-signature AUC comparison bar plot.")
    parser.add_argument("--config", required=True, help="Path to YAML config file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    input_csv = resolve_path(config["input_csv"], base_dir)
    output_dir = ensure_output_dir(resolve_path(config.get("output_dir", "figures/auc_barplot"), base_dir))

    dataset_col = config.get("dataset_col", "Dataset")
    signature_col = config.get("signature_col", "Signature")
    y_col = config.get("y_col", "AUC")
    auc_table = read_csv_checked(input_csv, required_columns=[dataset_col, signature_col, y_col])

    summary = plot_auc_barplot(
        auc_table=auc_table,
        output_dir=output_dir,
        signature_order=config.get("signature_order"),
        y_col=y_col,
        signature_col=signature_col,
        dataset_col=dataset_col,
        ylim=tuple(config.get("ylim", [0, 1.05])),
        order_by_mean_desc=bool(config.get("order_by_mean_desc", False)),
        auc_label_offset=float(config.get("auc_label_offset", 0.035)),
        dpi=int(config.get("dpi", 300)),
        basename=config.get("basename", "auc_barplot"),
    )
    print(summary.to_string(index=False))
    print(f"Saved AUC bar plot outputs to: {output_dir}")


if __name__ == "__main__":
    main()

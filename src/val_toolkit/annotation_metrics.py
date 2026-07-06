from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .labels import display_cell_type, display_method
from .validation import read_csv_checked, require_columns


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    annotation_csv: Path
    truth_col: str
    method_cols: list[str]


def _safe_divide(numerator: int | float, denominator: int | float) -> float:
    if denominator == 0:
        return np.nan
    return float(numerator) / float(denominator)


def compute_one_vs_rest_metrics(
    annotations: pd.DataFrame,
    truth_col: str,
    method_cols: Iterable[str],
    cell_types: Iterable[str],
    dataset_name: str,
    drop_missing_predictions: bool = True,
) -> pd.DataFrame:
    """
    Compute one-vs-rest precision, recall, F1, and accuracy for each method and cell type.

    Parameters
    ----------
    annotations:
        One row per cell. Must include a ground-truth label column and one prediction column per method.
    truth_col:
        Column containing the reference label.
    method_cols:
        Columns containing predicted labels from annotation methods.
    cell_types:
        Cell types to evaluate one-vs-rest.
    dataset_name:
        Dataset name to attach to the output.
    drop_missing_predictions:
        If True, rows with missing truth or missing prediction for the method are excluded.
    """
    method_cols = list(method_cols)
    cell_types = list(cell_types)
    require_columns(annotations, [truth_col, *method_cols], table_name="annotations")

    rows: list[dict[str, object]] = []

    for method_col in method_cols:
        if drop_missing_predictions:
            valid = annotations[truth_col].notna() & annotations[method_col].notna()
            current = annotations.loc[valid, [truth_col, method_col]].copy()
        else:
            current = annotations[[truth_col, method_col]].copy()
            current[truth_col] = current[truth_col].fillna("__missing_truth__")
            current[method_col] = current[method_col].fillna("__missing_prediction__")

        for cell_type in cell_types:
            y_true = current[truth_col].astype(str).eq(str(cell_type))
            y_pred = current[method_col].astype(str).eq(str(cell_type))

            tp = int((y_true & y_pred).sum())
            fp = int((~y_true & y_pred).sum())
            fn = int((y_true & ~y_pred).sum())
            tn = int((~y_true & ~y_pred).sum())

            precision = _safe_divide(tp, tp + fp)
            recall = _safe_divide(tp, tp + fn)
            f1 = _safe_divide(2 * precision * recall, precision + recall)
            accuracy = _safe_divide(tp + tn, tp + fp + fn + tn)

            metric_values = {
                "Precision": precision,
                "Recall": recall,
                "F1": f1,
                "Accuracy": accuracy,
            }
            for metric_name, metric_value in metric_values.items():
                rows.append(
                    {
                        "Dataset": dataset_name,
                        "Method": method_col,
                        "Method Display": display_method(method_col),
                        "Cell Type": cell_type,
                        "Cell Type Display": display_cell_type(cell_type),
                        "Metric": metric_name,
                        "Value": metric_value,
                        "TP": tp,
                        "FP": fp,
                        "FN": fn,
                        "TN": tn,
                        "N evaluated": int(len(current)),
                    }
                )

    return pd.DataFrame(rows)


def compute_metrics_from_dataset_configs(
    datasets: Iterable[DatasetConfig],
    cell_types: Iterable[str],
    drop_missing_predictions: bool = True,
) -> pd.DataFrame:
    """Compute benchmark metrics for multiple dataset configs and concatenate results."""
    all_metrics: list[pd.DataFrame] = []
    for dataset in datasets:
        required = [dataset.truth_col, *dataset.method_cols]
        annotations = read_csv_checked(dataset.annotation_csv, required_columns=required)
        metrics = compute_one_vs_rest_metrics(
            annotations=annotations,
            truth_col=dataset.truth_col,
            method_cols=dataset.method_cols,
            cell_types=cell_types,
            dataset_name=dataset.name,
            drop_missing_predictions=drop_missing_predictions,
        )
        all_metrics.append(metrics)

    if not all_metrics:
        raise ValueError("No datasets were provided.")

    return pd.concat(all_metrics, ignore_index=True)

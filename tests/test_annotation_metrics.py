from __future__ import annotations

import pandas as pd

from val_toolkit.annotation_metrics import compute_one_vs_rest_metrics


def test_compute_one_vs_rest_metrics_basic() -> None:
    annotations = pd.DataFrame(
        {
            "truth": ["B", "B", "CD4", "CD4"],
            "method": ["B", "CD4", "CD4", "CD4"],
        }
    )
    metrics = compute_one_vs_rest_metrics(
        annotations=annotations,
        truth_col="truth",
        method_cols=["method"],
        cell_types=["B", "CD4"],
        dataset_name="test",
    )

    f1_b = metrics.loc[
        (metrics["Metric"] == "F1") & (metrics["Cell Type"] == "B"),
        "Value",
    ].iloc[0]
    recall_cd4 = metrics.loc[
        (metrics["Metric"] == "Recall") & (metrics["Cell Type"] == "CD4"),
        "Value",
    ].iloc[0]

    assert round(f1_b, 6) == round(2 / 3, 6)
    assert recall_cd4 == 1.0

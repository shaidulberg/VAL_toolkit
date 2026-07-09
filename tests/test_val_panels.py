from __future__ import annotations

import numpy as np
import pandas as pd

from val_toolkit.val_panels import assign_ranked_bins, benjamini_hochberg, consensus_from_annotation_votes


def test_benjamini_hochberg_monotonic_examples():
    q = benjamini_hochberg([0.01, 0.04, 0.03, np.nan])
    assert np.isclose(q[0], 0.03)
    assert np.isclose(q[1], 0.04)
    assert np.isclose(q[2], 0.04)
    assert np.isnan(q[3])


def test_consensus_from_annotation_votes():
    table = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "m1": ["B cells", "NK"],
            "m2": ["B", "Monocytes"],
            "m3": ["CD4 T", "NK cell"],
            "c1": [0.9, 0.7],
            "c2": [0.8, 0.6],
            "c3": [0.2, 0.9],
        }
    )
    out = consensus_from_annotation_votes(
        table,
        annotation_columns=["m1", "m2", "m3"],
        confidence_columns={"m1": "c1", "m2": "c2", "m3": "c3"},
    )
    assert out.loc[0, "consensus_cell_type"] == "B"
    assert out.loc[0, "VAL_votes"] == 2
    assert out.loc[1, "consensus_cell_type"] == "NK"
    assert out.loc[1, "VAL_votes"] == 2


def test_assign_ranked_bins():
    table = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(10)],
            "consensus_cell_type": ["B"] * 10,
            "VAL_votes": [3] * 5 + [2] * 5,
            "mean_consensus_confidence": np.linspace(1, 0.1, 10),
            "VAL_rank_score": np.linspace(4, 1, 10),
        }
    )
    out = assign_ranked_bins(table, cell_types=["B"], n_bins=10)
    assert set(out["bin_index"]) == set(range(1, 11))
    assert out.iloc[0]["bin_label"] == "0-10"

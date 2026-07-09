from __future__ import annotations

import numpy as np
import pandas as pd

from val_toolkit.acs import AnnotationMethodConfig, add_acs_rank_columns, compute_acs_cell_table, parse_annotation_methods
from val_toolkit.annotation_harmonization import flatten_label_mapping, harmonize_label


def test_flatten_label_mapping_supports_grouped_yaml_style():
    mapping = {"B": ["B cell", "Memory B"], "CD4 T": ["CD4 T cell"]}
    flat = flatten_label_mapping(mapping)
    assert flat["B cell"] == "B"
    assert flat["Memory B"] == "B"
    assert flat["CD4 T cell"] == "CD4 T"


def test_harmonize_label_uses_default_synonyms_and_other():
    assert harmonize_label("Classical monocyte") == "Monocytes"
    assert harmonize_label("Natural killer cell") == "NK"
    assert harmonize_label("Unresolved stromal label") == "Other"


def test_parse_annotation_methods():
    methods = parse_annotation_methods(
        [
            {"name": "CellTypist", "label_column": "ct_label", "confidence_column": "ct_conf"},
            {"name": "SingleR", "label_column": "sr_label"},
        ]
    )
    assert methods[0].name == "CellTypist"
    assert methods[1].confidence_column is None


def test_compute_acs_cell_table_val_and_confidence():
    table = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "m1_label": ["B cell", "NK cell", "Unknown"],
            "m2_label": ["Memory B", "Monocyte", "Other"],
            "m3_label": ["CD4 T cell", "Natural killer cell", "Unclear"],
            "m1_conf": [0.9, 0.7, 0.2],
            "m2_conf": [0.8, 0.6, 0.3],
            "m3_conf": [0.2, 0.95, 0.4],
        }
    )
    methods = [
        AnnotationMethodConfig("m1", "m1_label", "m1_conf"),
        AnnotationMethodConfig("m2", "m2_label", "m2_conf"),
        AnnotationMethodConfig("m3", "m3_label", "m3_conf"),
    ]
    out = compute_acs_cell_table(table, methods=methods)
    assert out.loc[0, "consensus_cell_type"] == "B"
    assert out.loc[0, "VAL_votes"] == 2
    assert np.isclose(out.loc[0, "mean_consensus_confidence"], 0.85)
    assert out.loc[1, "consensus_cell_type"] == "NK"
    assert out.loc[1, "VAL_votes"] == 2
    assert out.loc[2, "consensus_cell_type"] == "Other"


def test_add_acs_rank_columns_sorts_by_val_then_confidence():
    table = pd.DataFrame(
        {
            "cell_id": ["low_conf", "val2", "val3", "mid_conf"],
            "consensus_cell_type": ["B", "B", "B", "B"],
            "VAL_votes": [3, 2, 3, 3],
            "mean_consensus_confidence": [0.2, 0.99, 0.9, 0.5],
            "VAL_rank_score": [3.2, 2.99, 3.9, 3.5],
        }
    )
    out = add_acs_rank_columns(table, major_cell_types=["B"])
    assert out.sort_values("acs_rank_within_cell_type")["cell_id"].tolist() == ["val3", "mid_conf", "low_conf", "val2"]

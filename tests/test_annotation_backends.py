from __future__ import annotations

import pandas as pd

from val_toolkit.annotation_backends import (
    BackendStatus,
    extract_celltypist_annotation_table,
    merge_annotation_tables,
    statuses_to_frame,
)


class FakeCellTypistResult:
    def __init__(self) -> None:
        self.predicted_labels = pd.DataFrame(
            {"predicted_labels": ["B cells", "CD4 T cells", "Unclear"]},
            index=["cell1", "cell2", "cell3"],
        )
        self.probability_matrix = pd.DataFrame(
            {
                "B cells": [0.91, 0.10, 0.20],
                "CD4 T cells": [0.05, 0.87, 0.30],
                "Monocytes": [0.04, 0.03, 0.50],
            },
            index=["cell1", "cell2", "cell3"],
        )


def test_extract_celltypist_annotation_table_uses_predicted_label_probability() -> None:
    table = extract_celltypist_annotation_table(
        FakeCellTypistResult(),
        cell_ids=["cell1", "cell2", "cell3"],
        prefix="celltypist",
        label_source="predicted_labels",
    )

    assert list(table.columns) == ["cell_id", "celltypist_label", "celltypist_confidence"]
    assert table.loc[0, "celltypist_label"] == "B cells"
    assert table.loc[0, "celltypist_confidence"] == 0.91
    assert table.loc[1, "celltypist_confidence"] == 0.87
    # The label "Unclear" is not in the probability matrix; fallback to row max.
    assert table.loc[2, "celltypist_confidence"] == 0.50


def test_merge_annotation_tables_outer_merges_by_cell_id() -> None:
    left = pd.DataFrame({"cell_id": ["c1", "c2"], "celltypist_label": ["B", "NK"]})
    right = pd.DataFrame({"cell_id": ["c2", "c3"], "singler_hpca_label": ["NK", "CD4 T"]})

    merged = merge_annotation_tables([left, right]).sort_values("cell_id").reset_index(drop=True)

    assert list(merged["cell_id"]) == ["c1", "c2", "c3"]
    assert pd.isna(merged.loc[0, "singler_hpca_label"])
    assert merged.loc[1, "celltypist_label"] == "NK"
    assert merged.loc[2, "singler_hpca_label"] == "CD4 T"


def test_statuses_to_frame_has_expected_columns() -> None:
    frame = statuses_to_frame(
        [
            BackendStatus(
                name="CellTypist",
                enabled=True,
                status="completed",
                message="ok",
                output_csv="celltypist.csv",
                label_column="celltypist_label",
                confidence_column="celltypist_confidence",
            )
        ]
    )

    assert frame.loc[0, "name"] == "CellTypist"
    assert frame.loc[0, "status"] == "completed"
    assert frame.loc[0, "label_column"] == "celltypist_label"

from val_toolkit.annotation_backends import (
    SingleRConfig,
    _sanitize_reference_name,
    run_singler_backend,
)


def test_singler_reference_name_sanitization() -> None:
    assert _sanitize_reference_name("HPCA") == "hpca"
    assert _sanitize_reference_name("BlueprintEncode") == "encode"
    assert _sanitize_reference_name("DatabaseImmuneCellExpression") == "dice"
    assert _sanitize_reference_name("Novershtern") == "hema"


def test_singler_backend_reports_missing_rscript() -> None:
    try:
        run_singler_backend(object(), config=SingleRConfig(rscript="definitely_missing_Rscript_for_test"))
    except Exception as exc:
        assert "Rscript" in str(exc) or "definitely_missing_Rscript_for_test" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected SingleR backend to fail when Rscript is missing")

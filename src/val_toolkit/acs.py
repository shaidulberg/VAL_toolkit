from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .annotation_harmonization import CANONICAL_ACS_CELL_TYPES, harmonize_label, validate_major_cell_types
from .validation import require_columns


@dataclass(frozen=True)
class AnnotationMethodConfig:
    """One annotation method contributing a vote to ACS/VAL."""

    name: str
    label_column: str
    confidence_column: str | None = None


def parse_annotation_methods(config_items: Sequence[Mapping[str, Any]]) -> list[AnnotationMethodConfig]:
    methods: list[AnnotationMethodConfig] = []
    for i, item in enumerate(config_items):
        if "name" not in item or "label_column" not in item:
            raise ValueError(
                "Each annotation_methods entry must contain at least 'name' and 'label_column'. "
                f"Problem entry index: {i}."
            )
        methods.append(
            AnnotationMethodConfig(
                name=str(item["name"]),
                label_column=str(item["label_column"]),
                confidence_column=(str(item["confidence_column"]) if item.get("confidence_column") else None),
            )
        )
    if not methods:
        raise ValueError("At least one annotation method is required for ACS annotation-column mode.")
    return methods


def _safe_method_name(name: str) -> str:
    return (
        str(name)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("-", "_")
        .replace(".", "_")
    )


def annotation_table_from_obs(adata: Any, methods: Sequence[AnnotationMethodConfig]) -> pd.DataFrame:
    """Extract per-method labels/confidences from adata.obs."""
    required = []
    for method in methods:
        required.append(method.label_column)
        if method.confidence_column:
            required.append(method.confidence_column)
    require_columns(adata.obs, required, table_name="adata.obs")

    obs = adata.obs.copy()
    obs = obs.reset_index(drop=True)
    obs.insert(0, "cell_id", adata.obs_names.astype(str))
    columns = ["cell_id"] + required
    return obs.loc[:, columns].copy()


def annotation_table_from_csv(
    annotation_csv: str | Path,
    methods: Sequence[AnnotationMethodConfig],
    cell_id_column: str = "cell_id",
) -> pd.DataFrame:
    """Load per-method labels/confidences from an external CSV."""
    path = Path(annotation_csv)
    if not path.exists():
        raise FileNotFoundError(f"Annotation CSV does not exist: {path}")
    table = pd.read_csv(path)
    required = [cell_id_column]
    for method in methods:
        required.append(method.label_column)
        if method.confidence_column:
            required.append(method.confidence_column)
    require_columns(table, required, table_name=str(path))
    out = table.loc[:, required].copy()
    if cell_id_column != "cell_id":
        out = out.rename(columns={cell_id_column: "cell_id"})
    out["cell_id"] = out["cell_id"].astype(str)
    return out


def compute_acs_cell_table(
    annotation_table: pd.DataFrame,
    methods: Sequence[AnnotationMethodConfig],
    label_mapping: Mapping[str, Any] | None = None,
    cell_id_column: str = "cell_id",
    major_cell_types: Sequence[str] = CANONICAL_ACS_CELL_TYPES,
) -> pd.DataFrame:
    """
    Compute manuscript-style ACS from multiple annotation-method calls.

    For each cell:
      1. Harmonize each method's raw label to a major immune lineage.
      2. Select the majority lineage as the consensus label.
      3. Define VAL as the number of methods assigning the consensus label.
      4. Define mean_consensus_confidence from only the supporting methods.
      5. Store ACS helper fields used for ranking cells by VAL descending and then by
         mean_consensus_confidence descending.
    """
    major_cell_types = validate_major_cell_types(major_cell_types)
    require_columns(annotation_table, [cell_id_column], table_name="annotation_table")
    for method in methods:
        required = [method.label_column]
        if method.confidence_column:
            required.append(method.confidence_column)
        require_columns(annotation_table, required, table_name="annotation_table")

    n_methods = len(methods)
    rows: list[dict[str, Any]] = []

    for _, source in annotation_table.iterrows():
        out_row: dict[str, Any] = {"cell_id": str(source[cell_id_column])}
        method_major_labels: list[str] = []
        method_confidences: list[float] = []

        for method in methods:
            safe_name = _safe_method_name(method.name)
            raw_label = source[method.label_column]
            major_label = harmonize_label(raw_label, label_mapping=label_mapping)
            if major_label not in set(major_cell_types) | {"Other"}:
                major_label = "Other"

            if method.confidence_column:
                confidence = pd.to_numeric(source[method.confidence_column], errors="coerce")
                confidence = float(confidence) if np.isfinite(confidence) else np.nan
            else:
                confidence = 1.0
            if not np.isfinite(confidence):
                confidence = 1.0

            method_major_labels.append(major_label)
            method_confidences.append(confidence)
            out_row[f"{safe_name}_raw_label"] = raw_label
            out_row[f"{safe_name}_major_label"] = major_label
            out_row[f"{safe_name}_confidence"] = confidence

        labels = np.asarray(method_major_labels, dtype=object)
        confidences = np.asarray(method_confidences, dtype=float)

        candidate_labels = [label for label in dict.fromkeys(labels.tolist()) if label != "Other"]
        if not candidate_labels:
            consensus = "Other"
            val_votes = int(np.sum(labels == "Other"))
        else:
            best_label = None
            best_key: tuple[int, float, int] = (-1, -np.inf, -1)
            for candidate in candidate_labels:
                support = labels == candidate
                # Tie-breaker: votes, then mean supporting confidence, then canonical order.
                canonical_order = major_cell_types.index(candidate) if candidate in major_cell_types else -1
                key = (int(support.sum()), float(np.nanmean(confidences[support])), -canonical_order)
                if key > best_key:
                    best_key = key
                    best_label = candidate
            consensus = str(best_label)
            val_votes = int(np.sum(labels == consensus))

        support_mask = labels == consensus
        mean_support_confidence = float(np.nanmean(confidences[support_mask])) if support_mask.any() else np.nan

        out_row.update(
            {
                "consensus_cell_type": consensus,
                "VAL": val_votes,
                "VAL_votes": val_votes,
                "VAL_fraction": val_votes / n_methods if n_methods else np.nan,
                "mean_consensus_confidence": mean_support_confidence,
                "ACS_rank_score": val_votes + mean_support_confidence / max(n_methods + 1, 1),
                # Backward-compatible alias used by the existing binning code.
                "VAL_rank_score": val_votes + mean_support_confidence / max(n_methods + 1, 1),
                "n_annotation_methods": n_methods,
                "n_supporting_methods": val_votes,
                "supporting_methods": ";".join(
                    method.name for method, label in zip(methods, labels, strict=False) if label == consensus
                ),
            }
        )
        rows.append(out_row)

    return pd.DataFrame(rows)


def add_acs_rank_columns(
    acs_cell_table: pd.DataFrame,
    major_cell_types: Sequence[str] = CANONICAL_ACS_CELL_TYPES,
) -> pd.DataFrame:
    """Add ACS rank/percentile columns within each consensus cell type."""
    require_columns(
        acs_cell_table,
        ["cell_id", "consensus_cell_type", "VAL_votes", "mean_consensus_confidence", "VAL_rank_score"],
        table_name="acs_cell_table",
    )
    pieces: list[pd.DataFrame] = []
    for cell_type in major_cell_types:
        current = acs_cell_table.loc[acs_cell_table["consensus_cell_type"].eq(cell_type)].copy()
        if current.empty:
            continue
        current = current.sort_values(
            ["VAL_votes", "mean_consensus_confidence", "VAL_rank_score", "cell_id"],
            ascending=[False, False, False, True],
        ).reset_index(drop=True)
        current["acs_rank_within_cell_type"] = np.arange(1, len(current) + 1)
        current["acs_percentile_within_cell_type"] = current["acs_rank_within_cell_type"] / len(current) * 100.0
        pieces.append(current)
    if not pieces:
        return acs_cell_table.copy()
    ranked = pd.concat(pieces, ignore_index=True)
    other = acs_cell_table.loc[~acs_cell_table["cell_id"].isin(ranked["cell_id"])].copy()
    if not other.empty:
        other["acs_rank_within_cell_type"] = np.nan
        other["acs_percentile_within_cell_type"] = np.nan
        ranked = pd.concat([ranked, other], ignore_index=True)
    return ranked

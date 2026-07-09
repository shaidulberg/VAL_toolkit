from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import mannwhitneyu, rankdata

from .labels import display_cell_type
from .validation import ensure_output_dir, require_columns

CANONICAL_CELL_TYPES: tuple[str, ...] = ("B", "CD4 T", "CD8 T", "DC", "Monocytes", "NK")

DEFAULT_MARKER_SETS: dict[str, tuple[str, ...]] = {
    "B": ("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "BANK1", "CD37", "HLA-DRA"),
    "CD4 T": ("CD3D", "CD3E", "CD3G", "TRAC", "CD4", "IL7R", "CCR7", "LTB"),
    "CD8 T": ("CD3D", "CD3E", "CD3G", "TRAC", "CD8A", "CD8B", "GZMK", "CCL5"),
    "DC": ("FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1", "IL3RA", "CLEC4C", "LILRA4", "LAMP3"),
    "Monocytes": ("LST1", "S100A8", "S100A9", "LYZ", "FCN1", "CTSS", "MS4A7", "FCGR3A"),
    "NK": ("NKG7", "GNLY", "PRF1", "KLRD1", "NCAM1", "GZMB", "FCGR3A", "KLRF1"),
}

LABEL_SYNONYMS: dict[str, str] = {
    "b": "B",
    "b cell": "B",
    "b cells": "B",
    "b_cell": "B",
    "b_cells": "B",
    "cd4": "CD4 T",
    "cd4 t": "CD4 T",
    "cd4 t cell": "CD4 T",
    "cd4 t cells": "CD4 T",
    "cd4_t": "CD4 T",
    "cd4_t_cells": "CD4 T",
    "cd8": "CD8 T",
    "cd8 t": "CD8 T",
    "cd8 t cell": "CD8 T",
    "cd8 t cells": "CD8 T",
    "cd8_t": "CD8 T",
    "cd8_t_cells": "CD8 T",
    "dc": "DC",
    "dendritic": "DC",
    "dendritic cell": "DC",
    "dendritic cells": "DC",
    "monocyte": "Monocytes",
    "monocytes": "Monocytes",
    "mono": "Monocytes",
    "myeloid/monocyte": "Monocytes",
    "nk": "NK",
    "nk cell": "NK",
    "nk cells": "NK",
    "natural killer": "NK",
    "natural killer cell": "NK",
    "natural killer cells": "NK",
}


@dataclass(frozen=True)
class SampleColumnConfig:
    patient_id: str
    response: str
    timepoint: str | None = None
    sample_id: str | None = None


def canonicalize_label(label: Any, user_label_map: Mapping[str, str] | None = None) -> str:
    """Map an arbitrary annotation label to one of the canonical immune labels or Other."""
    if pd.isna(label):
        return "Other"
    text = str(label).strip()
    if user_label_map and text in user_label_map:
        return str(user_label_map[text])
    lowered = text.lower().replace("-", " ").replace("_", " ").strip()
    if lowered in LABEL_SYNONYMS:
        return LABEL_SYNONYMS[lowered]
    for key, value in LABEL_SYNONYMS.items():
        if key in lowered:
            return value
    return "Other"


def _get_expression_matrix(adata: Any, layer: str | None = None, use_raw: bool = False) -> tuple[Any, pd.Index]:
    if use_raw:
        if getattr(adata, "raw", None) is None:
            raise ValueError("Config requested use_raw=true, but adata.raw is not present.")
        matrix = adata.raw.X
        var_names = pd.Index(adata.raw.var_names.astype(str))
    elif layer:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' was not found in adata.layers. Available layers: {list(adata.layers.keys())}")
        matrix = adata.layers[layer]
        var_names = pd.Index(adata.var_names.astype(str))
    else:
        matrix = adata.X
        var_names = pd.Index(adata.var_names.astype(str))
    return matrix, var_names


def _subset_matrix_columns(matrix: Any, column_indices: Sequence[int]) -> np.ndarray:
    sub = matrix[:, list(column_indices)]
    if sparse.issparse(sub):
        return sub.toarray()
    return np.asarray(sub)


def _maybe_log1p(values: np.ndarray, transform: str = "auto_log1p") -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if transform == "none":
        return values
    if transform == "log1p":
        return np.log1p(np.maximum(values, 0))
    if transform == "auto_log1p":
        finite = values[np.isfinite(values)]
        if finite.size and np.nanmax(finite) > 30:
            return np.log1p(np.maximum(values, 0))
        return values
    raise ValueError(f"Unknown expression transform: {transform}")


def _zscore_columns(values: np.ndarray) -> np.ndarray:
    means = np.nanmean(values, axis=0)
    stds = np.nanstd(values, axis=0)
    stds[stds == 0] = 1.0
    return (values - means) / stds


def _rank_percentile_columns(values: np.ndarray) -> np.ndarray:
    ranked = np.zeros_like(values, dtype=float)
    for j in range(values.shape[1]):
        ranked[:, j] = rankdata(values[:, j], method="average") / max(values.shape[0], 1)
    return ranked


def _softmax_confidence(score_frame: pd.DataFrame) -> pd.DataFrame:
    scores = score_frame.to_numpy(dtype=float)
    shifted = scores - np.nanmax(scores, axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    probs = exp_scores / np.nansum(exp_scores, axis=1, keepdims=True)
    return pd.DataFrame(probs, index=score_frame.index, columns=score_frame.columns)


def marker_ensemble_annotations(
    adata: Any,
    marker_sets: Mapping[str, Sequence[str]] | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    expression_transform: str = "auto_log1p",
    min_markers_per_cell_type: int = 2,
) -> pd.DataFrame:
    """
    Build three lightweight marker-based annotation votes directly from an AnnData object.

    This is the self-contained mode for users who only bring an h5ad file. It produces a
    VAL-like agreement score among three marker-scoring strategies. For manuscript-like VAL,
    use `consensus_from_obs_annotations` with labels from external annotation tools.
    """
    marker_sets = marker_sets or DEFAULT_MARKER_SETS
    matrix, var_names = _get_expression_matrix(adata, layer=layer, use_raw=use_raw)
    gene_to_idx = {gene.upper(): i for i, gene in enumerate(var_names.astype(str).str.upper())}
    cell_ids = pd.Index(adata.obs_names.astype(str), name="cell_id")

    mean_z_scores: dict[str, np.ndarray] = {}
    detection_scores: dict[str, np.ndarray] = {}
    percentile_scores: dict[str, np.ndarray] = {}
    markers_used: dict[str, list[str]] = {}

    for cell_type, genes in marker_sets.items():
        indices: list[int] = []
        used: list[str] = []
        for gene in genes:
            idx = gene_to_idx.get(str(gene).upper())
            if idx is not None:
                indices.append(idx)
                used.append(str(gene))
        markers_used[cell_type] = used
        if len(indices) < min_markers_per_cell_type:
            raise ValueError(
                f"Only {len(indices)} marker(s) found for '{cell_type}' but min_markers_per_cell_type="
                f"{min_markers_per_cell_type}. Found markers: {used}. Check gene symbols or lower the threshold."
            )
        expr = _subset_matrix_columns(matrix, indices)
        expr = _maybe_log1p(expr, transform=expression_transform)
        mean_z_scores[cell_type] = np.nanmean(_zscore_columns(expr), axis=1)
        detection_scores[cell_type] = np.nanmean(expr > 0, axis=1)
        percentile_scores[cell_type] = np.nanmean(_rank_percentile_columns(expr), axis=1)

    score_tables = {
        "marker_mean_z": pd.DataFrame(mean_z_scores, index=cell_ids),
        "marker_detection": pd.DataFrame(detection_scores, index=cell_ids),
        "marker_percentile": pd.DataFrame(percentile_scores, index=cell_ids),
    }

    rows = pd.DataFrame(index=cell_ids)
    for method_name, score_table in score_tables.items():
        conf_table = _softmax_confidence(score_table)
        labels = score_table.idxmax(axis=1)
        rows[f"{method_name}_label"] = labels.to_numpy()
        rows[f"{method_name}_confidence"] = conf_table.max(axis=1).to_numpy()
        rows[f"{method_name}_score"] = score_table.max(axis=1).to_numpy()

    rows.attrs["markers_used"] = markers_used
    return rows.reset_index()


def consensus_from_annotation_votes(
    annotation_table: pd.DataFrame,
    annotation_columns: Sequence[str],
    confidence_columns: Mapping[str, str] | None = None,
    label_map: Mapping[str, str] | None = None,
    cell_id_col: str = "cell_id",
) -> pd.DataFrame:
    """Compute consensus label, VAL votes, VAL fraction, and rank score from method labels."""
    confidence_columns = confidence_columns or {}
    require_columns(annotation_table, [cell_id_col, *annotation_columns], table_name="annotation_table")
    if confidence_columns:
        require_columns(annotation_table, confidence_columns.values(), table_name="annotation_table")

    rows: list[dict[str, Any]] = []
    n_methods = len(annotation_columns)
    if n_methods == 0:
        raise ValueError("At least one annotation column is required to compute VAL.")

    for _, row in annotation_table.iterrows():
        labels: list[str] = []
        confidences: list[float] = []
        per_method: dict[str, Any] = {cell_id_col: row[cell_id_col]}
        for method_col in annotation_columns:
            label = canonicalize_label(row[method_col], user_label_map=label_map)
            labels.append(label)
            conf_col = confidence_columns.get(method_col)
            confidence = float(pd.to_numeric(row[conf_col], errors="coerce")) if conf_col else 1.0
            if not np.isfinite(confidence):
                confidence = 1.0
            confidences.append(confidence)
            per_method[f"{method_col}_canonical"] = label
            per_method[f"{method_col}_confidence_used"] = confidence

        labels_array = np.asarray(labels, dtype=object)
        confidences_array = np.asarray(confidences, dtype=float)
        candidates = [label for label in dict.fromkeys(labels) if label != "Other"]
        if not candidates:
            consensus = "Other"
            val_votes = int(np.sum(labels_array == "Other"))
        else:
            best_label = None
            best_key = (-1, -np.inf)
            for candidate in candidates:
                mask = labels_array == candidate
                key = (int(mask.sum()), float(np.nanmean(confidences_array[mask])))
                if key > best_key:
                    best_key = key
                    best_label = candidate
            consensus = str(best_label)
            val_votes = int(best_key[0])

        support_mask = labels_array == consensus
        mean_support_conf = float(np.nanmean(confidences_array[support_mask])) if support_mask.any() else np.nan
        per_method.update(
            {
                "consensus_cell_type": consensus,
                "VAL_votes": val_votes,
                "VAL_fraction": val_votes / n_methods,
                "mean_consensus_confidence": mean_support_conf,
                "VAL_rank_score": val_votes + mean_support_conf / max(n_methods + 1, 1),
                "n_annotation_methods": n_methods,
            }
        )
        rows.append(per_method)

    return pd.DataFrame(rows)


def build_val_table_from_marker_ensemble(
    adata: Any,
    marker_sets: Mapping[str, Sequence[str]] | None = None,
    layer: str | None = None,
    use_raw: bool = False,
    expression_transform: str = "auto_log1p",
    min_markers_per_cell_type: int = 2,
) -> pd.DataFrame:
    annotations = marker_ensemble_annotations(
        adata=adata,
        marker_sets=marker_sets,
        layer=layer,
        use_raw=use_raw,
        expression_transform=expression_transform,
        min_markers_per_cell_type=min_markers_per_cell_type,
    )
    annotation_columns = ["marker_mean_z_label", "marker_detection_label", "marker_percentile_label"]
    confidence_columns = {
        "marker_mean_z_label": "marker_mean_z_confidence",
        "marker_detection_label": "marker_detection_confidence",
        "marker_percentile_label": "marker_percentile_confidence",
    }
    return consensus_from_annotation_votes(
        annotations,
        annotation_columns=annotation_columns,
        confidence_columns=confidence_columns,
        cell_id_col="cell_id",
    ).merge(annotations, on="cell_id", how="left")


def build_val_table_from_obs(
    adata: Any,
    annotation_columns: Sequence[str],
    confidence_columns: Mapping[str, str] | None = None,
    label_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    require_columns(adata.obs, annotation_columns, table_name="adata.obs")
    if confidence_columns:
        require_columns(adata.obs, confidence_columns.values(), table_name="adata.obs")
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", adata.obs_names.astype(str))
    return consensus_from_annotation_votes(
        obs,
        annotation_columns=annotation_columns,
        confidence_columns=confidence_columns,
        label_map=label_map,
        cell_id_col="cell_id",
    )


def assign_ranked_bins(
    val_table: pd.DataFrame,
    cell_types: Sequence[str] = CANONICAL_CELL_TYPES,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Assign non-overlapping ranked bins within each consensus cell type."""
    require_columns(
        val_table,
        ["cell_id", "consensus_cell_type", "VAL_votes", "mean_consensus_confidence", "VAL_rank_score"],
        table_name="val_table",
    )
    if n_bins <= 0:
        raise ValueError("n_bins must be positive.")

    pieces: list[pd.DataFrame] = []
    for cell_type in cell_types:
        current = val_table.loc[val_table["consensus_cell_type"].eq(cell_type)].copy()
        if current.empty:
            continue
        current = current.sort_values(
            ["VAL_votes", "mean_consensus_confidence", "VAL_rank_score", "cell_id"],
            ascending=[False, False, False, True],
        ).reset_index(drop=True)
        current["rank_within_cell_type"] = np.arange(1, len(current) + 1)
        current["rank_percentile"] = current["rank_within_cell_type"] / len(current) * 100.0
        bin_index = np.ceil(current["rank_percentile"] / (100.0 / n_bins)).astype(int)
        bin_index = bin_index.clip(lower=1, upper=n_bins)
        current["bin_index"] = bin_index
        current["bin_start_percent"] = (bin_index - 1) * (100.0 / n_bins)
        current["bin_end_percent"] = bin_index * (100.0 / n_bins)
        current["bin_label"] = current["bin_start_percent"].astype(int).astype(str) + "-" + current[
            "bin_end_percent"
        ].astype(int).astype(str)
        pieces.append(current)
    if not pieces:
        raise ValueError("No cells matched the requested cell_types after consensus labeling.")
    return pd.concat(pieces, ignore_index=True)


def _make_sample_table(adata: Any, sample_columns: SampleColumnConfig) -> pd.DataFrame:
    required = [sample_columns.patient_id, sample_columns.response]
    if sample_columns.timepoint:
        required.append(sample_columns.timepoint)
    if sample_columns.sample_id:
        required.append(sample_columns.sample_id)
    require_columns(adata.obs, required, table_name="adata.obs")

    obs = adata.obs[required].copy()
    obs.insert(0, "cell_id", adata.obs_names.astype(str))
    if sample_columns.sample_id:
        obs["sample_id"] = obs[sample_columns.sample_id].astype(str)
    elif sample_columns.timepoint:
        obs["sample_id"] = obs[sample_columns.patient_id].astype(str) + "__" + obs[sample_columns.timepoint].astype(str)
    else:
        obs["sample_id"] = obs[sample_columns.patient_id].astype(str)
    obs["patient_id"] = obs[sample_columns.patient_id].astype(str)
    obs["timepoint"] = obs[sample_columns.timepoint].astype(str) if sample_columns.timepoint else "NA"
    obs["response"] = obs[sample_columns.response].astype(str)
    return obs[["cell_id", "sample_id", "patient_id", "timepoint", "response"]].reset_index(drop=True)


def compute_ranked_bin_abundance(
    adata: Any,
    binned_cells: pd.DataFrame,
    sample_columns: SampleColumnConfig,
    cell_types: Sequence[str] = CANONICAL_CELL_TYPES,
    n_bins: int = 10,
    denominator: str = "all_cells",
) -> pd.DataFrame:
    """Compute one row per sample x cell-type x ranked bin with zero-filled proportions."""
    sample_by_cell = _make_sample_table(adata, sample_columns)
    sample_summary = sample_by_cell.groupby("sample_id", as_index=False).agg(
        patient_id=("patient_id", "first"),
        timepoint=("timepoint", "first"),
        response=("response", "first"),
        total_cells=("cell_id", "size"),
    )

    if denominator not in {"all_cells", "cell_type_cells"}:
        raise ValueError("denominator must be 'all_cells' or 'cell_type_cells'.")

    bins = pd.DataFrame(
        {
            "bin_index": np.arange(1, n_bins + 1),
            "bin_start_percent": (np.arange(1, n_bins + 1) - 1) * (100.0 / n_bins),
            "bin_end_percent": np.arange(1, n_bins + 1) * (100.0 / n_bins),
        }
    )
    bins["bin_label"] = bins["bin_start_percent"].astype(int).astype(str) + "-" + bins["bin_end_percent"].astype(
        int
    ).astype(str)
    grid = (
        sample_summary[["sample_id", "patient_id", "timepoint", "response", "total_cells"]]
        .assign(_key=1)
        .merge(pd.DataFrame({"consensus_cell_type": list(cell_types), "_key": 1}), on="_key")
        .merge(bins.assign(_key=1), on="_key")
        .drop(columns="_key")
    )

    binned_with_sample = binned_cells.merge(sample_by_cell[["cell_id", "sample_id"]], on="cell_id", how="left")
    counts = (
        binned_with_sample.groupby(["sample_id", "consensus_cell_type", "bin_index"], as_index=False)
        .size()
        .rename(columns={"size": "bin_cell_count"})
    )
    out = grid.merge(counts, on=["sample_id", "consensus_cell_type", "bin_index"], how="left")
    out["bin_cell_count"] = out["bin_cell_count"].fillna(0).astype(int)

    if denominator == "cell_type_cells":
        cell_type_totals = (
            binned_with_sample.groupby(["sample_id", "consensus_cell_type"], as_index=False)
            .size()
            .rename(columns={"size": "denominator_cells"})
        )
        out = out.merge(cell_type_totals, on=["sample_id", "consensus_cell_type"], how="left")
        out["denominator_cells"] = out["denominator_cells"].fillna(0).astype(int)
    else:
        out["denominator_cells"] = out["total_cells"].astype(int)

    out["proportion"] = np.where(out["denominator_cells"] > 0, out["bin_cell_count"] / out["denominator_cells"], np.nan)
    return out.sort_values(["consensus_cell_type", "bin_index", "sample_id"]).reset_index(drop=True)


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    p_valid = p[valid]
    order = np.argsort(p_valid)
    ranked = p_valid[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out_valid = np.empty_like(adjusted)
    out_valid[order] = adjusted
    q[valid] = out_valid
    return q


def compute_ranked_bin_associations(
    abundance: pd.DataFrame,
    positive_response: str,
    negative_response: str,
    fdr_scope: str = "global",
) -> pd.DataFrame:
    require_columns(
        abundance,
        ["consensus_cell_type", "bin_index", "bin_end_percent", "bin_label", "sample_id", "response", "proportion"],
        table_name="abundance",
    )
    rows: list[dict[str, Any]] = []
    for (cell_type, bin_index), current in abundance.groupby(["consensus_cell_type", "bin_index"], sort=True):
        pos = pd.to_numeric(current.loc[current["response"].eq(str(positive_response)), "proportion"], errors="coerce").dropna()
        neg = pd.to_numeric(current.loc[current["response"].eq(str(negative_response)), "proportion"], errors="coerce").dropna()
        if len(pos) == 0 or len(neg) == 0:
            p_value = np.nan
            statistic = np.nan
        elif np.nanstd(pd.concat([pos, neg])) == 0:
            p_value = 1.0
            statistic = np.nan
        else:
            result = mannwhitneyu(pos, neg, alternative="two-sided")
            p_value = float(result.pvalue)
            statistic = float(result.statistic)
        meta = current.iloc[0]
        rows.append(
            {
                "consensus_cell_type": cell_type,
                "bin_index": int(bin_index),
                "bin_start_percent": float(meta["bin_start_percent"]),
                "bin_end_percent": float(meta["bin_end_percent"]),
                "bin_label": meta["bin_label"],
                "positive_response": positive_response,
                "negative_response": negative_response,
                "n_positive_samples": int(len(pos)),
                "n_negative_samples": int(len(neg)),
                "mean_positive_proportion": float(pos.mean()) if len(pos) else np.nan,
                "mean_negative_proportion": float(neg.mean()) if len(neg) else np.nan,
                "effect_positive_minus_negative": (float(pos.mean()) - float(neg.mean())) if len(pos) and len(neg) else np.nan,
                "mannwhitneyu_statistic": statistic,
                "p_value": p_value,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    if fdr_scope == "global":
        out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    elif fdr_scope == "within_cell_type":
        out["q_value"] = np.nan
        for cell_type, idx in out.groupby("consensus_cell_type").groups.items():
            out.loc[idx, "q_value"] = benjamini_hochberg(out.loc[idx, "p_value"].to_numpy())
    else:
        raise ValueError("fdr_scope must be 'global' or 'within_cell_type'.")
    positive = out["q_value"].replace(0, np.nan).dropna()
    zero_replacement = positive.min() * 0.1 if not positive.empty else 1e-300
    q_for_plot = out["q_value"].replace(0, zero_replacement)
    out["neg_log10_q_value"] = -np.log10(q_for_plot)
    out["neg_log10_q_value"] = out["neg_log10_q_value"].replace([np.inf, -np.inf], np.nan).fillna(0)
    return out.sort_values(["consensus_cell_type", "bin_index"]).reset_index(drop=True)


def _save_figure(fig: plt.Figure, output_dir: Path, basename: str, dpi: int = 300) -> None:
    output_dir = ensure_output_dir(output_dir)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(output_dir / f"{basename}.{extension}", dpi=dpi, bbox_inches="tight")


def plot_ranked_bin_panels(
    association_results: pd.DataFrame,
    output_dir: str | Path,
    cell_types: Sequence[str] = CANONICAL_CELL_TYPES,
    fdr_threshold: float = 0.05,
    dpi: int = 300,
    basename: str = "ranked_bin_panels_A_to_F",
    x_axis_prefix: str = "Ranked",
) -> None:
    require_columns(
        association_results,
        ["consensus_cell_type", "bin_end_percent", "q_value", "neg_log10_q_value"],
        table_name="association_results",
    )
    output_dir = ensure_output_dir(output_dir)
    y_threshold = -np.log10(fdr_threshold)
    max_y = max(float(association_results["neg_log10_q_value"].max()), y_threshold)
    ylim_top = max(2.2, max_y + 0.35)

    fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.6), constrained_layout=True)
    axes_flat = axes.ravel()
    letters = list("ABCDEF")

    for ax, letter, cell_type in zip(axes_flat, letters, cell_types):
        current = association_results.loc[association_results["consensus_cell_type"].eq(cell_type)].copy()
        current = current.sort_values("bin_end_percent")
        if current.empty:
            ax.text(0.5, 0.5, f"No {display_cell_type(cell_type)} cells", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        x = current["bin_end_percent"].to_numpy(dtype=float)
        y = current["neg_log10_q_value"].to_numpy(dtype=float)
        scatter = ax.scatter(x, y, c=y, s=58, cmap="viridis", edgecolor="black", linewidth=0.25, zorder=3)
        ax.axhline(y_threshold, color="red", linestyle=":", linewidth=1.0, alpha=0.8)
        best = current.sort_values(["q_value", "bin_end_percent"], ascending=[True, True]).iloc[0]
        best_x = float(best["bin_end_percent"])
        ax.axvline(best_x, color="black", linestyle="--", linewidth=0.8, alpha=0.65)
        ax.text(best_x, ylim_top * 0.94, f"{int(round(best_x))}%", ha="center", va="top", fontsize=10)
        ax.set_xlim(0, 105)
        ax.set_ylim(0, ylim_top)
        ax.set_xticks(np.arange(0, 101, 10))
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.35)
        ax.set_xlabel(f"{x_axis_prefix} {display_cell_type(cell_type)} cell bin (%)", fontsize=11)
        ax.set_ylabel("-log10(FDR q-value)", fontsize=11)
        ax.set_title(letter, loc="left", fontweight="bold", fontsize=15)
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
        cbar.ax.set_ylabel("-log10(FDR q-value)", fontsize=7)
        cbar.ax.tick_params(labelsize=7)

    _save_figure(fig, output_dir, basename=basename, dpi=dpi)
    plt.close(fig)

    # Individual panels are useful for users who want to assemble their own figure.
    for letter, cell_type in zip(letters, cell_types):
        fig_i, ax_i = plt.subplots(figsize=(4.6, 3.6))
        current = association_results.loc[association_results["consensus_cell_type"].eq(cell_type)].copy()
        current = current.sort_values("bin_end_percent")
        if current.empty:
            ax_i.text(0.5, 0.5, f"No {display_cell_type(cell_type)} cells", ha="center", va="center", transform=ax_i.transAxes)
            ax_i.set_axis_off()
        else:
            x = current["bin_end_percent"].to_numpy(dtype=float)
            y = current["neg_log10_q_value"].to_numpy(dtype=float)
            scatter = ax_i.scatter(x, y, c=y, s=62, cmap="viridis", edgecolor="black", linewidth=0.25, zorder=3)
            ax_i.axhline(y_threshold, color="red", linestyle=":", linewidth=1.0, alpha=0.8)
            best = current.sort_values(["q_value", "bin_end_percent"], ascending=[True, True]).iloc[0]
            best_x = float(best["bin_end_percent"])
            ax_i.axvline(best_x, color="black", linestyle="--", linewidth=0.8, alpha=0.65)
            ax_i.text(best_x, ylim_top * 0.94, f"{int(round(best_x))}%", ha="center", va="top", fontsize=10)
            ax_i.set_xlim(0, 105)
            ax_i.set_ylim(0, ylim_top)
            ax_i.set_xticks(np.arange(0, 101, 10))
            ax_i.grid(True, linestyle="--", linewidth=0.4, alpha=0.35)
            ax_i.set_xlabel(f"{x_axis_prefix} {display_cell_type(cell_type)} cell bin (%)", fontsize=11)
            ax_i.set_ylabel("-log10(FDR q-value)", fontsize=11)
            ax_i.set_title(letter, loc="left", fontweight="bold", fontsize=15)
            cbar = fig_i.colorbar(scatter, ax=ax_i, fraction=0.046, pad=0.03)
            cbar.ax.set_ylabel("-log10(FDR q-value)", fontsize=7)
            cbar.ax.tick_params(labelsize=7)
        _save_figure(fig_i, output_dir, basename=f"ranked_bin_panel_{letter}_{cell_type.replace(' ', '_')}", dpi=dpi)
        plt.close(fig_i)

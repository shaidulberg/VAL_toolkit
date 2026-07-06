from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import sem

from .labels import display_cell_type, display_method
from .validation import ensure_output_dir, require_columns


def _mean_sem(values: pd.Series) -> tuple[float, float]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if values.empty:
        return np.nan, np.nan
    if len(values) == 1:
        return float(values.iloc[0]), 0.0
    return float(values.mean()), float(sem(values, nan_policy="omit"))


def _jitter_positions(center: float, n: int, width: float = 0.16) -> np.ndarray:
    if n <= 1:
        return np.array([center])
    return np.linspace(center - width / 2, center + width / 2, n)


def _save_figure(fig: plt.Figure, output_dir: Path, basename: str, dpi: int = 300) -> None:
    output_dir = ensure_output_dir(output_dir)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(output_dir / f"{basename}.{extension}", dpi=dpi, bbox_inches="tight")


def plot_annotation_metric_by_method(
    metrics: pd.DataFrame,
    output_dir: str | Path,
    metric: str = "F1",
    cell_type_order: Sequence[str] | None = None,
    method_order: Sequence[str] | None = None,
    ylim: tuple[float, float] = (0, 1.05),
    dpi: int = 300,
) -> None:
    """Create one plot per metric: x-axis methods, grouped by cell type."""
    require_columns(metrics, ["Metric", "Value", "Method", "Cell Type", "Dataset"], table_name="metrics")
    output_dir = ensure_output_dir(output_dir)

    df = metrics.loc[metrics["Metric"].eq(metric)].copy()
    if df.empty:
        raise ValueError(f"No rows found for metric: {metric}")

    df["Method Display"] = df["Method"].map(display_method)
    df["Cell Type Display"] = df["Cell Type"].map(display_cell_type)

    if method_order is None:
        method_order = list(dict.fromkeys(df["Method"].tolist()))
    if cell_type_order is None:
        cell_type_order = list(dict.fromkeys(df["Cell Type"].tolist()))

    x = np.arange(len(method_order), dtype=float)
    n_groups = len(cell_type_order)
    bar_width = min(0.78 / max(n_groups, 1), 0.12)

    fig, ax = plt.subplots(figsize=(max(10, len(method_order) * 0.95), 6.2))

    for i, cell_type in enumerate(cell_type_order):
        offset = (i - (n_groups - 1) / 2) * bar_width
        means: list[float] = []
        errors: list[float] = []
        for method in method_order:
            values = df.loc[(df["Method"].eq(method)) & (df["Cell Type"].eq(cell_type)), "Value"]
            mean_value, sem_value = _mean_sem(values)
            means.append(mean_value)
            errors.append(sem_value)

        positions = x + offset
        ax.bar(
            positions,
            means,
            width=bar_width,
            yerr=errors,
            capsize=2.5,
            label=display_cell_type(cell_type),
            edgecolor="black",
            linewidth=0.6,
        )

        for method_index, method in enumerate(method_order):
            values = df.loc[(df["Method"].eq(method)) & (df["Cell Type"].eq(cell_type)), "Value"].dropna()
            jitter = _jitter_positions(positions[method_index], len(values))
            ax.scatter(jitter, values, s=18, edgecolor="black", linewidth=0.35, zorder=3)

    ax.set_ylabel(metric, fontsize=14)
    ax.set_ylim(*ylim)
    ax.set_xticks(x)
    ax.set_xticklabels([display_method(method) for method in method_order], rotation=45, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.legend(frameon=False, fontsize=10, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    _save_figure(fig, output_dir, f"annotation_benchmark_by_method_{metric}", dpi=dpi)
    plt.close(fig)


def plot_annotation_metric_by_cell_type(
    metrics: pd.DataFrame,
    output_dir: str | Path,
    metric: str = "F1",
    cell_type_order: Sequence[str] | None = None,
    method_order: Sequence[str] | None = None,
    ylim: tuple[float, float] = (0, 1.05),
    dpi: int = 300,
) -> None:
    """Create one plot per metric: x-axis cell types, grouped by method."""
    require_columns(metrics, ["Metric", "Value", "Method", "Cell Type", "Dataset"], table_name="metrics")
    output_dir = ensure_output_dir(output_dir)

    df = metrics.loc[metrics["Metric"].eq(metric)].copy()
    if df.empty:
        raise ValueError(f"No rows found for metric: {metric}")

    if method_order is None:
        method_order = list(dict.fromkeys(df["Method"].tolist()))
    if cell_type_order is None:
        cell_type_order = list(dict.fromkeys(df["Cell Type"].tolist()))

    x = np.arange(len(cell_type_order), dtype=float)
    n_groups = len(method_order)
    bar_width = min(0.78 / max(n_groups, 1), 0.12)

    fig, ax = plt.subplots(figsize=(max(10, len(cell_type_order) * 1.35), 6.2))

    for i, method in enumerate(method_order):
        offset = (i - (n_groups - 1) / 2) * bar_width
        means: list[float] = []
        errors: list[float] = []
        for cell_type in cell_type_order:
            values = df.loc[(df["Method"].eq(method)) & (df["Cell Type"].eq(cell_type)), "Value"]
            mean_value, sem_value = _mean_sem(values)
            means.append(mean_value)
            errors.append(sem_value)

        positions = x + offset
        ax.bar(
            positions,
            means,
            width=bar_width,
            yerr=errors,
            capsize=2.5,
            label=display_method(method),
            edgecolor="black",
            linewidth=0.6,
        )

        for cell_index, cell_type in enumerate(cell_type_order):
            values = df.loc[(df["Method"].eq(method)) & (df["Cell Type"].eq(cell_type)), "Value"].dropna()
            jitter = _jitter_positions(positions[cell_index], len(values))
            ax.scatter(jitter, values, s=18, edgecolor="black", linewidth=0.35, zorder=3)

    ax.set_ylabel(metric, fontsize=14)
    ax.set_ylim(*ylim)
    ax.set_xticks(x)
    ax.set_xticklabels([display_cell_type(cell_type) for cell_type in cell_type_order], rotation=45, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.legend(frameon=False, fontsize=9, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    _save_figure(fig, output_dir, f"annotation_benchmark_by_cell_type_{metric}", dpi=dpi)
    plt.close(fig)


def plot_auc_barplot(
    auc_table: pd.DataFrame,
    output_dir: str | Path,
    signature_order: Sequence[str] | None = None,
    y_col: str = "AUC",
    signature_col: str = "Signature",
    dataset_col: str = "Dataset",
    ylim: tuple[float, float] = (0, 1.05),
    order_by_mean_desc: bool = False,
    auc_label_offset: float = 0.035,
    dpi: int = 300,
    basename: str = "auc_barplot",
) -> pd.DataFrame:
    """Plot mean ± SEM AUC by signature with individual dataset points."""
    require_columns(auc_table, [dataset_col, signature_col, y_col], table_name="auc_table")
    output_dir = ensure_output_dir(output_dir)

    df = auc_table[[dataset_col, signature_col, y_col]].copy()
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col, signature_col, dataset_col])

    if signature_order is None:
        if order_by_mean_desc:
            signature_order = (
                df.groupby(signature_col)[y_col]
                .mean()
                .sort_values(ascending=False)
                .index.tolist()
            )
        else:
            signature_order = list(dict.fromkeys(df[signature_col].tolist()))

    summary_rows: list[dict[str, object]] = []
    for signature in signature_order:
        values = df.loc[df[signature_col].eq(signature), y_col]
        mean_value, sem_value = _mean_sem(values)
        summary_rows.append(
            {
                signature_col: signature,
                "Mean AUC": mean_value,
                "SEM": sem_value,
                "N datasets": int(values.notna().sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)

    x = np.arange(len(signature_order), dtype=float)
    fig, ax = plt.subplots(figsize=(max(10, len(signature_order) * 0.9), 6.2))

    ax.bar(
        x,
        summary["Mean AUC"].to_numpy(dtype=float),
        yerr=summary["SEM"].to_numpy(dtype=float),
        capsize=3,
        edgecolor="black",
        linewidth=0.7,
    )

    for i, signature in enumerate(signature_order):
        values = df.loc[df[signature_col].eq(signature), y_col].dropna().reset_index(drop=True)
        jitter = _jitter_positions(float(i), len(values), width=0.18)
        ax.scatter(jitter, values, s=28, edgecolor="black", linewidth=0.45, zorder=3)

        mean_value = summary.loc[summary[signature_col].eq(signature), "Mean AUC"].iloc[0]
        if pd.notna(mean_value):
            ax.text(
                i,
                min(mean_value + auc_label_offset, ylim[1] - 0.015),
                f"{mean_value:.3f}",
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold",
            )

    ax.set_ylabel("AUC", fontsize=14)
    ax.set_ylim(*ylim)
    ax.set_xticks(x)
    ax.set_xticklabels(signature_order, rotation=45, ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    _save_figure(fig, output_dir, basename, dpi=dpi)
    summary.to_csv(output_dir / f"{basename}_summary.csv", index=False)
    df.to_csv(output_dir / f"{basename}_raw.csv", index=False)
    plt.close(fig)
    return summary

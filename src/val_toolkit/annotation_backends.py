from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

try:
    from scipy import sparse
except Exception:  # pragma: no cover - scipy is a package dependency, this is defensive.
    sparse = None  # type: ignore[assignment]


@dataclass(frozen=True)
class BackendStatus:
    """Execution status for one optional annotation backend."""

    name: str
    enabled: bool
    status: str
    message: str
    output_csv: str | None = None
    label_column: str | None = None
    confidence_column: str | None = None


class BackendUnavailableError(RuntimeError):
    """Raised when an optional annotation backend cannot be run in this environment."""


@dataclass(frozen=True)
class CellTypistConfig:
    """Configuration for the optional CellTypist annotation backend."""

    model: str = "Immune_All_Low.pkl"
    majority_voting: bool = False
    over_clustering: str | None = None
    use_gpu: bool = False
    prefix: str = "celltypist"
    label_source: str = "predicted_labels"
    download_models: bool = False
    force_update_models: bool = False


@dataclass(frozen=True)
class SingleRConfig:
    """
    Configuration placeholder for a future automatic SingleR backend.

    The current toolkit does not automatically execute SingleR. It records a clear
    skip/status message and points users to the annotation-column ACS workflow.
    """

    references: tuple[str, ...] = ("HPCA", "Monaco", "DICE", "BlueprintEncode", "DatabaseImmuneCellExpression", "Novershtern")
    prefix: str = "singler"
    rscript: str = "Rscript"
    on_unavailable: str = "skip"


@dataclass(frozen=True)
class SCimilarityConfig:
    """
    Configuration placeholder for a future automatic SCimilarity backend.

    The current toolkit does not automatically execute SCimilarity because users
    must provide a compatible SCimilarity installation and model path.
    """

    model_path: str | None = None
    prefix: str = "scimilarity"
    on_unavailable: str = "skip"


def _as_1d_numeric(values: Any) -> np.ndarray:
    """Convert dense or sparse row/column summaries to a flat numeric array."""
    if sparse is not None and sparse.issparse(values):
        values = values.A
    arr = np.asarray(values).reshape(-1)
    return arr.astype(float, copy=False)


def _copy_anndata(adata: Any) -> Any:
    try:
        return adata.copy()
    except Exception as exc:  # pragma: no cover
        raise TypeError("Expected an AnnData-like object with a .copy() method.") from exc


def prepare_celltypist_input(
    adata: Any,
    *,
    layer: str | None = None,
    use_raw: bool = False,
    normalize: bool = False,
    target_sum: float = 10000.0,
    log1p: bool = False,
) -> Any:
    """
    Return an AnnData copy prepared for CellTypist.

    CellTypist expects a cell-by-gene matrix with gene symbols in ``var_names``.
    If ``normalize`` and ``log1p`` are False, this function leaves expression values
    unchanged. This is safest when the user already stores log1p-normalized expression
    in ``adata.X`` or the selected layer.
    """
    out = _copy_anndata(adata)

    if use_raw:
        if out.raw is None:
            raise ValueError("use_raw=True was requested, but adata.raw is None.")
        out = out.raw.to_adata()
    elif layer:
        if layer not in out.layers:
            raise ValueError(f"Requested CellTypist layer '{layer}' is not present in adata.layers.")
        out.X = out.layers[layer].copy()

    if normalize:
        sums = _as_1d_numeric(out.X.sum(axis=1))
        scale = np.divide(target_sum, sums, out=np.zeros_like(sums, dtype=float), where=sums > 0)
        if sparse is not None and sparse.issparse(out.X):
            out.X = sparse.diags(scale).dot(out.X)
        else:
            out.X = np.asarray(out.X, dtype=float) * scale[:, None]

    if log1p:
        if sparse is not None and sparse.issparse(out.X):
            out.X = out.X.copy()
            out.X.data = np.log1p(out.X.data)
        else:
            out.X = np.log1p(np.asarray(out.X, dtype=float))

    return out


def load_celltypist() -> Any:
    """Import CellTypist lazily so the core ACS workflow works without it installed."""
    try:
        import celltypist  # type: ignore
    except ImportError as exc:
        raise BackendUnavailableError(
            "CellTypist is not installed in this environment. Install the optional backend with:\n"
            "  python -m pip install celltypist\n"
            "or recreate the environment after adding CellTypist to the dependencies."
        ) from exc
    return celltypist


def maybe_download_celltypist_models(download_models: bool, force_update: bool) -> None:
    """Optionally download CellTypist built-in models."""
    if not download_models:
        return
    try:
        from celltypist import models  # type: ignore
    except ImportError as exc:
        raise BackendUnavailableError("CellTypist model download requested, but CellTypist is not installed.") from exc
    models.download_models(force_update=force_update)


def _choose_celltypist_label_column(labels: pd.DataFrame, label_source: str, majority_voting: bool) -> str:
    if label_source == "auto":
        if majority_voting and "majority_voting" in labels.columns:
            return "majority_voting"
        return "predicted_labels" if "predicted_labels" in labels.columns else str(labels.columns[0])
    if label_source not in labels.columns:
        raise ValueError(
            f"Requested CellTypist label_source '{label_source}' was not found in predicted_labels columns: "
            f"{list(labels.columns)}"
        )
    return label_source


def extract_celltypist_annotation_table(
    result: Any,
    *,
    cell_ids: list[str] | pd.Index,
    prefix: str = "celltypist",
    label_source: str = "auto",
    majority_voting: bool = False,
) -> pd.DataFrame:
    """
    Convert a CellTypist AnnotationResult to a compact annotation table.

    The confidence assigned to each cell is the CellTypist probability for the
    predicted label whenever that label exists in the probability matrix. If the
    selected label is unavailable in the probability matrix, the row maximum
    probability is used as a conservative fallback.
    """
    labels = result.predicted_labels.copy()
    probabilities = result.probability_matrix.copy()

    cell_ids = pd.Index([str(x) for x in cell_ids], name="cell_id")
    labels.index = labels.index.astype(str)
    probabilities.index = probabilities.index.astype(str)
    probabilities.columns = probabilities.columns.astype(str)

    chosen_col = _choose_celltypist_label_column(labels, label_source=label_source, majority_voting=majority_voting)
    chosen_labels = labels[chosen_col].astype(str).reindex(cell_ids)
    prob = probabilities.reindex(cell_ids)

    max_prob = prob.max(axis=1).astype(float)
    confidences: list[float] = []
    for cell_id, label in chosen_labels.items():
        if pd.isna(label):
            confidences.append(np.nan)
        elif str(label) in prob.columns:
            value = prob.loc[cell_id, str(label)]
            confidences.append(float(value) if pd.notna(value) else np.nan)
        else:
            value = max_prob.loc[cell_id]
            confidences.append(float(value) if pd.notna(value) else np.nan)

    return pd.DataFrame(
        {
            "cell_id": cell_ids.astype(str),
            f"{prefix}_label": chosen_labels.to_numpy(dtype=object),
            f"{prefix}_confidence": confidences,
        }
    )


def run_celltypist_backend(
    adata: Any,
    *,
    config: CellTypistConfig | None = None,
    expression_layer: str | None = None,
    use_raw: bool = False,
    normalize: bool = False,
    target_sum: float = 10000.0,
    log1p: bool = False,
) -> tuple[pd.DataFrame, Any]:
    """Run CellTypist and return a compact annotation table plus the raw result."""
    cfg = config or CellTypistConfig()
    maybe_download_celltypist_models(cfg.download_models, cfg.force_update_models)
    celltypist = load_celltypist()

    query = prepare_celltypist_input(
        adata,
        layer=expression_layer,
        use_raw=use_raw,
        normalize=normalize,
        target_sum=target_sum,
        log1p=log1p,
    )

    result = celltypist.annotate(
        query,
        model=cfg.model,
        majority_voting=cfg.majority_voting,
        over_clustering=cfg.over_clustering,
        use_GPU=cfg.use_gpu,
    )
    table = extract_celltypist_annotation_table(
        result,
        cell_ids=list(query.obs_names.astype(str)),
        prefix=cfg.prefix,
        label_source=cfg.label_source,
        majority_voting=cfg.majority_voting,
    )
    return table, result


def check_singler_available(rscript: str = "Rscript") -> tuple[bool, str]:
    """Return whether Rscript appears available for a future SingleR bridge."""
    if shutil.which(rscript) is None:
        return False, f"{rscript!r} was not found on PATH."
    try:
        probe = subprocess.run(
            [rscript, "-e", "cat(requireNamespace('SingleR', quietly=TRUE))"],
            check=False,
            text=True,
            capture_output=True,
            timeout=30,
        )
    except Exception as exc:  # pragma: no cover - depends on host R setup.
        return False, f"Could not probe SingleR through {rscript!r}: {exc}"
    if "TRUE" not in probe.stdout:
        return False, "R is available, but the Bioconductor package SingleR is not installed."
    return True, "Rscript and SingleR are available. Automatic SingleR execution is scaffolded for a future release."


def run_singler_backend(*_: Any, config: SingleRConfig | None = None, **__: Any) -> tuple[pd.DataFrame, Any]:
    """
    Placeholder for automatic SingleR execution.

    Use ``scripts/run_acs_figures.py`` with annotation columns or an annotation CSV
    to include SingleR outputs today. This function intentionally raises a clear
    message rather than silently producing incomplete manuscript-style VAL results.
    """
    cfg = config or SingleRConfig()
    available, message = check_singler_available(cfg.rscript)
    if not available:
        raise BackendUnavailableError(message)
    raise NotImplementedError(
        "Automatic SingleR execution is not implemented yet. Export SingleR labels/scores "
        "into adata.obs or an annotation CSV and include them in configs/acs_figures.example.yaml."
    )


def check_scimilarity_available(model_path: str | None = None) -> tuple[bool, str]:
    """Return whether the SCimilarity Python package and optional model path appear available."""
    try:
        import scimilarity  # noqa: F401  # type: ignore
    except ImportError:
        return False, "The Python package scimilarity is not installed."
    if model_path and not Path(model_path).expanduser().exists():
        return False, f"SCimilarity model_path does not exist: {model_path}"
    return True, "SCimilarity package is importable. Automatic SCimilarity execution is scaffolded for a future release."


def run_scimilarity_backend(*_: Any, config: SCimilarityConfig | None = None, **__: Any) -> tuple[pd.DataFrame, Any]:
    """
    Placeholder for automatic SCimilarity execution.

    Use the annotation-column ACS workflow with precomputed SCimilarity outputs today.
    """
    cfg = config or SCimilarityConfig()
    available, message = check_scimilarity_available(cfg.model_path)
    if not available:
        raise BackendUnavailableError(message)
    raise NotImplementedError(
        "Automatic SCimilarity execution is not implemented yet. Export SCimilarity labels/confidences "
        "into adata.obs or an annotation CSV and include them in configs/acs_figures.example.yaml."
    )


def attach_annotation_table_to_obs(adata: Any, table: pd.DataFrame) -> Any:
    """Attach annotation label/confidence columns to an AnnData object's obs."""
    if "cell_id" not in table.columns:
        raise ValueError("annotation table must contain a 'cell_id' column")
    out = _copy_anndata(adata)
    keyed = table.copy()
    keyed["cell_id"] = keyed["cell_id"].astype(str)
    keyed = keyed.set_index("cell_id")
    keyed = keyed.reindex(out.obs_names.astype(str))
    for col in keyed.columns:
        out.obs[col] = keyed[col].to_numpy()
    return out


def merge_annotation_tables(tables: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Outer-merge compact backend annotation tables by cell_id."""
    tables = [t.copy() for t in tables if t is not None and not t.empty]
    if not tables:
        return pd.DataFrame(columns=["cell_id"])
    merged = tables[0]
    if "cell_id" not in merged.columns:
        raise ValueError("Every annotation table must contain a 'cell_id' column.")
    for table in tables[1:]:
        if "cell_id" not in table.columns:
            raise ValueError("Every annotation table must contain a 'cell_id' column.")
        overlap = (set(merged.columns) & set(table.columns)) - {"cell_id"}
        if overlap:
            raise ValueError(f"Annotation table columns overlap unexpectedly: {sorted(overlap)}")
        merged = merged.merge(table, on="cell_id", how="outer")
    return merged


def statuses_to_frame(statuses: Iterable[BackendStatus]) -> pd.DataFrame:
    """Convert backend status objects to a CSV-ready table."""
    return pd.DataFrame([s.__dict__ for s in statuses])


def parse_backend_enabled(config: Mapping[str, Any], name: str, default: bool = False) -> bool:
    """Read backend enabled status from a config mapping."""
    section = config.get(name, {})
    if isinstance(section, Mapping):
        return bool(section.get("enabled", default))
    return bool(section)

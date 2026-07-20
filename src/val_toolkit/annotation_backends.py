from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

try:
    from scipy import sparse
    from scipy.io import mmwrite
except Exception:  # pragma: no cover - scipy is a package dependency, this is defensive.
    sparse = None  # type: ignore[assignment]
    mmwrite = None  # type: ignore[assignment]


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
    """Configuration for the optional automatic SingleR backend.

    SingleR is executed through an Rscript bridge because SingleR and celldex are
    Bioconductor packages. The backend exports the selected AnnData expression
    matrix to MatrixMarket format, runs SingleR for each requested celldex
    reference, and imports compact label/confidence columns for ACS.
    """

    references: tuple[str, ...] = ("HPCA", "Monaco", "DICE", "BlueprintEncode", "Novershtern")
    prefix: str = "singler"
    rscript: str = "Rscript"
    label_field: str = "label.main"
    prune_labels: bool = True
    score_column: str = "assigned_label_score"
    assay_type_ref: str = "logcounts"
    keep_temp: bool = False
    on_unavailable: str = "skip"


@dataclass(frozen=True)
class SCimilarityConfig:
    """Configuration for the optional automatic SCimilarity backend.

    SCimilarity is executed through its Python API. Users must install the
    optional ``scimilarity`` package and provide a compatible local model
    directory. The backend aligns the query AnnData object to the model gene
    order, computes embeddings, predicts nearest-neighbor cell-type labels, and
    exports compact label/confidence columns for ACS.
    """

    model_path: str | None = None
    prefix: str = "scimilarity"
    use_gpu: bool = False
    k: int = 50
    ef: int = 100
    weighting: bool = False
    confidence_column: str = "vsAll"
    gene_overlap_threshold: int = 5000
    safelist: tuple[str, ...] = ()
    blocklist: tuple[str, ...] = ()
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
    """Return whether Rscript, SingleR, celldex, and Matrix are available."""
    if shutil.which(rscript) is None:
        return False, f"{rscript!r} was not found on PATH."
    probe_expr = "cat(all(sapply(c('SingleR','celldex','Matrix','SummarizedExperiment'), requireNamespace, quietly=TRUE)))"
    try:
        probe = subprocess.run(
            [rscript, "-e", probe_expr],
            check=False,
            text=True,
            capture_output=True,
            timeout=30,
        )
    except Exception as exc:  # pragma: no cover - depends on host R setup.
        return False, f"Could not probe SingleR dependencies through {rscript!r}: {exc}"
    if "TRUE" not in probe.stdout:
        details = (probe.stderr or probe.stdout or "").strip()
        return (
            False,
            "R is available, but one or more required Bioconductor packages are missing: "
            "SingleR, celldex, Matrix, SummarizedExperiment. Install them with BiocManager."
            + (f" Probe output: {details}" if details else ""),
        )
    return True, "Rscript and required SingleR dependencies are available."


def _sanitize_reference_name(reference: str) -> str:
    """Convert a reference name to a stable lowercase suffix for output columns."""
    cleaned = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(reference)).strip("_")
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    aliases = {
        "humanprimarycellatlas": "hpca",
        "databaseimmunecellexpression": "dice",
        "databaseimmunecellexpressiondata": "dice",
        "blueprintencode": "encode",
        "blueprint": "encode",
        "novershternhematopoietic": "hema",
        "novershtern": "hema",
    }
    return aliases.get(cleaned, cleaned or "reference")


def _write_lines(path: Path, values: Iterable[Any]) -> None:
    with path.open("w") as handle:
        for value in values:
            handle.write(str(value) + "\n")


def export_anndata_matrix_for_r(
    adata: Any,
    output_dir: Path,
    *,
    expression_layer: str | None = None,
    use_raw: bool = False,
    normalize: bool = False,
    target_sum: float = 10000.0,
    log1p: bool = False,
) -> dict[str, Path]:
    """Export an AnnData expression matrix for the R SingleR bridge.

    The MatrixMarket file is written as cells x genes. The R bridge transposes it
    to genes x cells before calling SingleR.
    """
    if mmwrite is None or sparse is None:
        raise BackendUnavailableError("scipy.io.mmwrite is required to export matrices for the SingleR backend.")

    output_dir.mkdir(parents=True, exist_ok=True)
    query = prepare_celltypist_input(
        adata,
        layer=expression_layer,
        use_raw=use_raw,
        normalize=normalize,
        target_sum=target_sum,
        log1p=log1p,
    )

    matrix_path = output_dir / "query_cells_by_genes.mtx"
    cell_ids_path = output_dir / "cell_ids.txt"
    gene_ids_path = output_dir / "gene_ids.txt"

    x = query.X
    if sparse.issparse(x):
        matrix_to_write = x.tocoo()
    else:
        matrix_to_write = sparse.coo_matrix(np.asarray(x, dtype=float))
    mmwrite(str(matrix_path), matrix_to_write)
    _write_lines(cell_ids_path, query.obs_names.astype(str))
    _write_lines(gene_ids_path, query.var_names.astype(str))

    return {"matrix": matrix_path, "cell_ids": cell_ids_path, "gene_ids": gene_ids_path}


def _singler_r_script_path() -> Path:
    return Path(__file__).resolve().parents[2] / "scripts" / "run_singler_backend.R"


def _run_single_reference_singler(
    *,
    cfg: SingleRConfig,
    exported: dict[str, Path],
    reference: str,
    output_csv: Path,
) -> pd.DataFrame:
    script_path = _singler_r_script_path()
    if not script_path.exists():
        raise FileNotFoundError(f"SingleR R bridge script not found: {script_path}")

    cmd = [
        cfg.rscript,
        str(script_path),
        "--matrix",
        str(exported["matrix"]),
        "--cell_ids",
        str(exported["cell_ids"]),
        "--gene_ids",
        str(exported["gene_ids"]),
        "--output_csv",
        str(output_csv),
        "--reference",
        str(reference),
        "--label_field",
        cfg.label_field,
        "--prune_labels",
        str(cfg.prune_labels).lower(),
        "--score_column",
        cfg.score_column,
        "--assay_type_ref",
        cfg.assay_type_ref,
    ]
    result = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "SingleR R bridge failed for reference "
            f"{reference!r}.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    table = pd.read_csv(output_csv)
    if not {"cell_id", "label", "confidence"}.issubset(table.columns):
        raise ValueError(f"SingleR output for {reference!r} is missing required columns: {output_csv}")
    suffix = _sanitize_reference_name(reference)
    out = pd.DataFrame(
        {
            "cell_id": table["cell_id"].astype(str),
            f"{cfg.prefix}_{suffix}_label": table["label"].astype(str),
            f"{cfg.prefix}_{suffix}_confidence": pd.to_numeric(table["confidence"], errors="coerce"),
        }
    )
    if "delta_next" in table.columns:
        out[f"{cfg.prefix}_{suffix}_delta_next"] = pd.to_numeric(table["delta_next"], errors="coerce")
    return out


def run_singler_backend(
    adata: Any,
    *,
    config: SingleRConfig | None = None,
    expression_layer: str | None = None,
    use_raw: bool = False,
    normalize: bool = False,
    target_sum: float = 10000.0,
    log1p: bool = False,
    work_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run SingleR through an Rscript bridge and return compact annotation columns."""
    cfg = config or SingleRConfig()
    available, message = check_singler_available(cfg.rscript)
    if not available:
        raise BackendUnavailableError(message)

    temp_ctx = None
    if work_dir is None:
        temp_ctx = tempfile.TemporaryDirectory(prefix="val_singler_")
        run_dir = Path(temp_ctx.name)
    else:
        run_dir = Path(work_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

    try:
        exported = export_anndata_matrix_for_r(
            adata,
            run_dir,
            expression_layer=expression_layer,
            use_raw=use_raw,
            normalize=normalize,
            target_sum=target_sum,
            log1p=log1p,
        )
        reference_tables: list[pd.DataFrame] = []
        reference_outputs: dict[str, str] = {}
        for reference in cfg.references:
            suffix = _sanitize_reference_name(reference)
            output_csv = run_dir / f"singler_{suffix}_raw.csv"
            table = _run_single_reference_singler(
                cfg=cfg,
                exported=exported,
                reference=reference,
                output_csv=output_csv,
            )
            reference_tables.append(table)
            reference_outputs[str(reference)] = str(output_csv)

        merged = merge_annotation_tables(reference_tables)
        raw = {"reference_outputs": reference_outputs, "work_dir": str(run_dir)}
        return merged, raw
    finally:
        if temp_ctx is not None and not cfg.keep_temp:
            temp_ctx.cleanup()


def check_scimilarity_available(model_path: str | None = None) -> tuple[bool, str]:
    """Return whether the SCimilarity Python package and model path appear available."""
    try:
        import scimilarity  # noqa: F401  # type: ignore
        from scimilarity.cell_annotation import CellAnnotation  # noqa: F401  # type: ignore
        from scimilarity.utils import align_dataset  # noqa: F401  # type: ignore
    except ImportError:
        return (
            False,
            "The Python package scimilarity is not installed or its annotation API is unavailable. "
            "Install it with `python -m pip install scimilarity` or from the Genentech/scimilarity source.",
        )
    if not model_path:
        return False, "SCimilarity requires a local model_path pointing to downloaded SCimilarity model files."
    if not Path(model_path).expanduser().exists():
        return False, f"SCimilarity model_path does not exist: {model_path}"
    return True, "SCimilarity package and model path are available."


def _load_scimilarity_api() -> tuple[Any, Any]:
    """Import SCimilarity annotation objects lazily."""
    try:
        from scimilarity.cell_annotation import CellAnnotation  # type: ignore
        from scimilarity.utils import align_dataset  # type: ignore
    except ImportError as exc:
        raise BackendUnavailableError(
            "The Python package scimilarity is not installed or does not expose the expected "
            "CellAnnotation/align_dataset API. Install it with `python -m pip install scimilarity` "
            "or from the Genentech/scimilarity source."
        ) from exc
    return CellAnnotation, align_dataset


def _choose_scimilarity_confidence(stats: pd.DataFrame, requested: str, weighting: bool) -> pd.Series:
    """Choose a numeric confidence column from SCimilarity KNN statistics."""
    preferred = requested
    if requested == "auto":
        preferred = "vsAll_weighted" if weighting else "vsAll"
    candidates = [preferred, "vsAll_weighted", "vsAll", "vs2nd_weighted", "vs2nd"]
    for col in candidates:
        if col in stats.columns:
            return pd.to_numeric(stats[col], errors="coerce")
    if "min_dist" in stats.columns:
        dist = pd.to_numeric(stats["min_dist"], errors="coerce")
        return 1.0 / (1.0 + dist)
    return pd.Series(np.nan, index=stats.index, dtype=float)


def extract_scimilarity_annotation_table(
    *,
    predictions: Any,
    stats: pd.DataFrame,
    cell_ids: list[str] | pd.Index,
    prefix: str = "scimilarity",
    confidence_column: str = "vsAll",
    weighting: bool = False,
) -> pd.DataFrame:
    """Convert SCimilarity KNN predictions/statistics into compact ACS columns."""
    cell_ids = pd.Index([str(x) for x in cell_ids], name="cell_id")
    if isinstance(predictions, pd.Series):
        labels = predictions.copy()
        labels.index = labels.index.astype(str)
        labels = labels.reindex(cell_ids)
    else:
        labels = pd.Series(np.asarray(predictions).reshape(-1), index=cell_ids)

    stats = stats.copy()
    if len(stats.index) == len(cell_ids):
        stats.index = cell_ids
    else:
        stats.index = stats.index.astype(str)
        stats = stats.reindex(cell_ids)

    confidence = _choose_scimilarity_confidence(stats, confidence_column, weighting).reindex(cell_ids)
    out = pd.DataFrame(
        {
            "cell_id": cell_ids.astype(str),
            f"{prefix}_label": labels.astype(object).to_numpy(),
            f"{prefix}_confidence": confidence.to_numpy(dtype=float),
        }
    )
    for optional_col in ["min_dist", "max_dist", "vs2nd", "vsAll", "vs2nd_weighted", "vsAll_weighted"]:
        if optional_col in stats.columns:
            out[f"{prefix}_{optional_col}"] = pd.to_numeric(stats[optional_col], errors="coerce").to_numpy()
    return out


def run_scimilarity_backend(
    adata: Any,
    *,
    config: SCimilarityConfig | None = None,
    expression_layer: str | None = None,
    use_raw: bool = False,
    normalize: bool = False,
    target_sum: float = 10000.0,
    log1p: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run SCimilarity cell-type annotation and return compact ACS columns."""
    cfg = config or SCimilarityConfig()
    available, message = check_scimilarity_available(cfg.model_path)
    if not available:
        raise BackendUnavailableError(message)

    CellAnnotation, align_dataset = _load_scimilarity_api()
    model_path = str(Path(str(cfg.model_path)).expanduser())
    annotator = CellAnnotation(model_path=model_path, use_gpu=cfg.use_gpu)

    if cfg.safelist:
        annotator.safelist_celltypes(list(cfg.safelist))
    if cfg.blocklist:
        annotator.blocklist_celltypes(list(cfg.blocklist))

    query = prepare_celltypist_input(
        adata,
        layer=expression_layer,
        use_raw=use_raw,
        normalize=normalize,
        target_sum=target_sum,
        log1p=log1p,
    )
    aligned = align_dataset(
        query,
        annotator.gene_order,
        gene_overlap_threshold=int(cfg.gene_overlap_threshold),
    )
    embeddings = annotator.get_embeddings(aligned.X)
    predictions, nn_idxs, nn_dists, stats = annotator.get_predictions_knn(
        embeddings,
        k=int(cfg.k),
        ef=int(cfg.ef),
        weighting=bool(cfg.weighting),
    )
    table = extract_scimilarity_annotation_table(
        predictions=predictions,
        stats=pd.DataFrame(stats),
        cell_ids=list(aligned.obs_names.astype(str)),
        prefix=cfg.prefix,
        confidence_column=cfg.confidence_column,
        weighting=cfg.weighting,
    )
    raw = {
        "model_path": model_path,
        "n_cells": int(aligned.n_obs),
        "n_genes_aligned": int(aligned.n_vars),
        "k": int(cfg.k),
        "ef": int(cfg.ef),
        "nn_idxs_shape": tuple(np.asarray(nn_idxs).shape),
        "nn_dists_shape": tuple(np.asarray(nn_dists).shape),
    }
    return table, raw


def _coerce_label_value_for_obs(value: Any) -> str | float:
    """Return a plain string or np.nan for AnnData obs label storage.

    Some optional backends, especially SCimilarity, may return labels as
    numpy arrays, tuples, or other Python objects. anndata/h5py cannot write
    those objects directly as variable-length strings, so label columns are
    normalized before being attached to ``adata.obs``.
    """
    if value is None:
        return np.nan
    if isinstance(value, float) and np.isnan(value):
        return np.nan
    if isinstance(value, (list, tuple, set, np.ndarray)):
        arr = np.asarray(list(value) if isinstance(value, set) else value, dtype=object).reshape(-1)
        cleaned = [_coerce_label_value_for_obs(v) for v in arr]
        cleaned = [str(v) for v in cleaned if not (isinstance(v, float) and np.isnan(v)) and str(v) != ""]
        if not cleaned:
            return np.nan
        # Most prediction APIs return a single label. If multiple labels are
        # returned, keep a deterministic semicolon-separated representation.
        return ";".join(dict.fromkeys(cleaned))
    try:
        if bool(pd.isna(value)):
            return np.nan
    except Exception:
        pass
    return str(value)


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
        series = keyed[col]
        if col.endswith("_label"):
            out.obs[col] = pd.Categorical(series.map(_coerce_label_value_for_obs))
        else:
            out.obs[col] = series.to_numpy()
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

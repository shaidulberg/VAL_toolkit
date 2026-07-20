#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from val_toolkit.annotation_backends import (
    BackendStatus,
    BackendUnavailableError,
    CellTypistConfig,
    SCimilarityConfig,
    SingleRConfig,
    attach_annotation_table_to_obs,
    merge_annotation_tables,
    run_celltypist_backend,
    run_scimilarity_backend,
    run_singler_backend,
    statuses_to_frame,
)
from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.validation import ensure_output_dir


def _load_anndata():
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "This workflow requires anndata. Install it or recreate the conda environment from environment.yml."
        ) from exc
    return ad


def _coalesce(cli_value: Any, config_value: Any, default: Any = None) -> Any:
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return default


def _should_skip(section: dict[str, Any], exc: Exception) -> bool:
    on_unavailable = str(section.get("on_unavailable", "skip")).lower()
    if isinstance(exc, NotImplementedError):
        return on_unavailable != "error"
    if isinstance(exc, BackendUnavailableError):
        return on_unavailable != "error"
    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run optional annotation backends for downstream ACS figures. "
            "CellTypist, SingleR, and SCimilarity are supported as optional backends."
        )
    )
    parser.add_argument("--config", required=True, help="YAML config file.")
    parser.add_argument("--h5ad", default=None, help="Optional override for input_h5ad in the config.")
    parser.add_argument("--output-dir", default=None, help="Optional override for output_dir in the config.")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    input_h5ad = resolve_path(_coalesce(args.h5ad, config.get("input_h5ad")), base_dir=base_dir)
    output_dir = ensure_output_dir(
        resolve_path(_coalesce(args.output_dir, config.get("output_dir"), "figures/annotation_backends"), base_dir=base_dir)
    )

    ad = _load_anndata()
    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad does not exist: {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad)
    print(f"Loaded h5ad: {input_h5ad}")
    print(f"AnnData shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    expression_cfg = config.get("expression", {})
    backend_cfg = config.get("annotation_backends", {})
    output_cfg = config.get("outputs", {})

    tables = []
    statuses: list[BackendStatus] = []

    celltypist_cfg = backend_cfg.get("celltypist", {})
    if bool(celltypist_cfg.get("enabled", False)):
        print("\n[CellTypist] Running CellTypist backend...")
        prefix = str(celltypist_cfg.get("prefix", "celltypist"))
        try:
            table, _ = run_celltypist_backend(
                adata,
                config=CellTypistConfig(
                    model=str(celltypist_cfg.get("model", "Immune_All_Low.pkl")),
                    majority_voting=bool(celltypist_cfg.get("majority_voting", False)),
                    over_clustering=celltypist_cfg.get("over_clustering"),
                    use_gpu=bool(celltypist_cfg.get("use_gpu", False)),
                    prefix=prefix,
                    label_source=str(celltypist_cfg.get("label_source", "auto")),
                    download_models=bool(celltypist_cfg.get("download_models", False)),
                    force_update_models=bool(celltypist_cfg.get("force_update_models", False)),
                ),
                expression_layer=expression_cfg.get("layer"),
                use_raw=bool(expression_cfg.get("use_raw", False)),
                normalize=bool(expression_cfg.get("normalize_to_target_sum", False)),
                target_sum=float(expression_cfg.get("target_sum", 10000.0)),
                log1p=bool(expression_cfg.get("log1p", False)),
            )
            table_path = output_dir / f"{prefix}_annotations.csv"
            table.to_csv(table_path, index=False)
            tables.append(table)
            statuses.append(
                BackendStatus(
                    name="CellTypist",
                    enabled=True,
                    status="completed",
                    message="CellTypist completed successfully.",
                    output_csv=str(table_path),
                    label_column=f"{prefix}_label",
                    confidence_column=f"{prefix}_confidence",
                )
            )
            print(f"[CellTypist] Saved: {table_path}")
        except Exception as exc:
            if _should_skip(celltypist_cfg, exc):
                statuses.append(BackendStatus("CellTypist", True, "skipped", str(exc)))
                print(f"[CellTypist] Skipped: {exc}")
            else:
                raise
    else:
        statuses.append(BackendStatus("CellTypist", False, "disabled", "CellTypist backend disabled in config."))

    singler_cfg = backend_cfg.get("singler", {})
    if bool(singler_cfg.get("enabled", False)):
        print("\n[SingleR] Running SingleR backend through Rscript...")
        prefix = str(singler_cfg.get("prefix", "singler"))
        references = tuple(singler_cfg.get("references", ["HPCA", "Monaco", "DICE", "BlueprintEncode", "Novershtern"]))
        try:
            table, raw = run_singler_backend(
                adata,
                config=SingleRConfig(
                    references=references,
                    prefix=prefix,
                    rscript=str(singler_cfg.get("rscript", "Rscript")),
                    label_field=str(singler_cfg.get("label_field", "label.main")),
                    prune_labels=bool(singler_cfg.get("prune_labels", True)),
                    score_column=str(singler_cfg.get("score_column", "assigned_label_score")),
                    assay_type_ref=str(singler_cfg.get("assay_type_ref", "logcounts")),
                    keep_temp=bool(singler_cfg.get("keep_temp", False)),
                    on_unavailable=str(singler_cfg.get("on_unavailable", "skip")),
                ),
                expression_layer=expression_cfg.get("layer"),
                use_raw=bool(expression_cfg.get("use_raw", False)),
                normalize=bool(expression_cfg.get("normalize_to_target_sum", False)),
                target_sum=float(expression_cfg.get("target_sum", 10000.0)),
                log1p=bool(expression_cfg.get("log1p", False)),
                work_dir=output_dir / "singler_work" if bool(singler_cfg.get("keep_temp", False)) else None,
            )
            table_path = output_dir / f"{prefix}_annotations.csv"
            table.to_csv(table_path, index=False)
            tables.append(table)
            label_cols = [c for c in table.columns if c.endswith("_label")]
            confidence_cols = [c for c in table.columns if c.endswith("_confidence")]
            statuses.append(
                BackendStatus(
                    name="SingleR",
                    enabled=True,
                    status="completed",
                    message=f"SingleR completed successfully for references: {', '.join(references)}.",
                    output_csv=str(table_path),
                    label_column=",".join(label_cols),
                    confidence_column=",".join(confidence_cols),
                )
            )
            print(f"[SingleR] Saved: {table_path}")
        except Exception as exc:
            if _should_skip(singler_cfg, exc):
                statuses.append(BackendStatus("SingleR", True, "skipped", str(exc)))
                print(f"[SingleR] Skipped: {exc}")
            else:
                raise
    else:
        statuses.append(BackendStatus("SingleR", False, "disabled", "SingleR backend disabled in config."))

    scimilarity_cfg = backend_cfg.get("scimilarity", {})
    if bool(scimilarity_cfg.get("enabled", False)):
        print("\n[SCimilarity] Running SCimilarity backend...")
        prefix = str(scimilarity_cfg.get("prefix", "scimilarity"))
        try:
            table, raw = run_scimilarity_backend(
                adata,
                config=SCimilarityConfig(
                    model_path=scimilarity_cfg.get("model_path"),
                    prefix=prefix,
                    use_gpu=bool(scimilarity_cfg.get("use_gpu", False)),
                    k=int(scimilarity_cfg.get("k", 50)),
                    ef=int(scimilarity_cfg.get("ef", 100)),
                    weighting=bool(scimilarity_cfg.get("weighting", False)),
                    confidence_column=str(scimilarity_cfg.get("confidence_column", "vsAll")),
                    gene_overlap_threshold=int(scimilarity_cfg.get("gene_overlap_threshold", 5000)),
                    safelist=tuple(scimilarity_cfg.get("safelist", []) or []),
                    blocklist=tuple(scimilarity_cfg.get("blocklist", []) or []),
                    on_unavailable=str(scimilarity_cfg.get("on_unavailable", "skip")),
                ),
                expression_layer=expression_cfg.get("layer"),
                use_raw=bool(expression_cfg.get("use_raw", False)),
                normalize=bool(expression_cfg.get("normalize_to_target_sum", False)),
                target_sum=float(expression_cfg.get("target_sum", 10000.0)),
                log1p=bool(expression_cfg.get("log1p", False)),
            )
            table_path = output_dir / f"{prefix}_annotations.csv"
            table.to_csv(table_path, index=False)
            tables.append(table)
            statuses.append(
                BackendStatus(
                    name="SCimilarity",
                    enabled=True,
                    status="completed",
                    message=f"SCimilarity completed successfully with k={scimilarity_cfg.get('k', 50)}.",
                    output_csv=str(table_path),
                    label_column=f"{prefix}_label",
                    confidence_column=f"{prefix}_confidence",
                )
            )
            print(f"[SCimilarity] Saved: {table_path}")
        except Exception as exc:
            if _should_skip(scimilarity_cfg, exc):
                statuses.append(BackendStatus("SCimilarity", True, "skipped", str(exc)))
                print(f"[SCimilarity] Skipped: {exc}")
            else:
                raise
    else:
        statuses.append(BackendStatus("SCimilarity", False, "disabled", "SCimilarity backend disabled in config."))

    combined = merge_annotation_tables(tables)
    combined_csv = output_dir / output_cfg.get("combined_annotation_csv", "annotation_backend_columns.csv")
    combined.to_csv(combined_csv, index=False)
    print(f"\nSaved combined annotation table: {combined_csv}")

    status_table = statuses_to_frame(statuses)
    status_csv = output_dir / output_cfg.get("status_csv", "annotation_backend_status.csv")
    status_table.to_csv(status_csv, index=False)
    print(f"Saved backend status table: {status_csv}")

    if bool(output_cfg.get("write_h5ad", True)):
        annotated = attach_annotation_table_to_obs(adata, combined) if not combined.empty else adata.copy()
        h5ad_path = output_dir / output_cfg.get("annotated_h5ad", "annotation_backends_annotated.h5ad")
        annotated.write_h5ad(h5ad_path)
        print(f"Saved h5ad with completed backend columns in obs: {h5ad_path}")

    print("\nBackend status summary:")
    print(status_table.to_string(index=False))
    print(
        "\nNext step: include completed label/confidence columns in configs/acs_figures.example.yaml "
        "and run scripts/run_acs_figures.py to generate ACS figures A-F."
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from val_toolkit.annotation_backends import CellTypistConfig, attach_annotation_table_to_obs, run_celltypist_backend
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


def _coalesce(cli_value, config_value, default=None):
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return default


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the optional CellTypist backend and save annotation columns for downstream ACS figures."
    )
    parser.add_argument("--config", required=True, help="YAML config file.")
    parser.add_argument("--h5ad", default=None, help="Optional override for input_h5ad in the config.")
    parser.add_argument("--output-dir", default=None, help="Optional override for output_dir in the config.")
    parser.add_argument("--model", default=None, help="Optional override for celltypist.model.")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    input_h5ad = resolve_path(_coalesce(args.h5ad, config.get("input_h5ad")), base_dir=base_dir)
    output_dir = ensure_output_dir(
        resolve_path(_coalesce(args.output_dir, config.get("output_dir"), "figures/celltypist_annotations"), base_dir=base_dir)
    )

    celltypist_cfg = config.get("celltypist", {})
    expression_cfg = config.get("expression", {})
    output_cfg = config.get("outputs", {})

    ad = _load_anndata()
    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad does not exist: {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad)
    print(f"Loaded h5ad: {input_h5ad}")
    print(f"AnnData shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    cfg = CellTypistConfig(
        model=str(_coalesce(args.model, celltypist_cfg.get("model"), "Immune_All_Low.pkl")),
        majority_voting=bool(celltypist_cfg.get("majority_voting", False)),
        over_clustering=celltypist_cfg.get("over_clustering"),
        use_gpu=bool(celltypist_cfg.get("use_gpu", False)),
        prefix=str(celltypist_cfg.get("prefix", "celltypist")),
        label_source=str(celltypist_cfg.get("label_source", "auto")),
        download_models=bool(celltypist_cfg.get("download_models", False)),
        force_update_models=bool(celltypist_cfg.get("force_update_models", False)),
    )

    table, _ = run_celltypist_backend(
        adata,
        config=cfg,
        expression_layer=expression_cfg.get("layer"),
        use_raw=bool(expression_cfg.get("use_raw", False)),
        normalize=bool(expression_cfg.get("normalize_to_target_sum", False)),
        target_sum=float(expression_cfg.get("target_sum", 10000.0)),
        log1p=bool(expression_cfg.get("log1p", False)),
    )

    table_path = output_dir / output_cfg.get("annotation_csv", "celltypist_annotations.csv")
    table.to_csv(table_path, index=False)
    print(f"Saved CellTypist annotation table: {table_path}")

    if bool(output_cfg.get("write_h5ad", True)):
        annotated = attach_annotation_table_to_obs(adata, table)
        h5ad_path = output_dir / output_cfg.get("annotated_h5ad", "celltypist_annotated.h5ad")
        annotated.write_h5ad(h5ad_path)
        print(f"Saved h5ad with CellTypist columns in obs: {h5ad_path}")

    print("\nCellTypist columns added:")
    for col in table.columns:
        if col != "cell_id":
            print(f"  - {col}")


if __name__ == "__main__":
    main()

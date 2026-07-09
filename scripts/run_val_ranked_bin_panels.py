#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.val_panels import (
    CANONICAL_CELL_TYPES,
    SampleColumnConfig,
    assign_ranked_bins,
    build_val_table_from_marker_ensemble,
    build_val_table_from_obs,
    compute_ranked_bin_abundance,
    compute_ranked_bin_associations,
    plot_ranked_bin_panels,
)
from val_toolkit.validation import ensure_output_dir


def _load_anndata() -> Any:
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "This workflow requires anndata. Install the package with `pip install anndata` "
            "or recreate the conda environment from environment.yml."
        ) from exc
    return ad


def _coalesce(cli_value: str | None, config_value: Any, default: Any = None) -> Any:
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return default


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate VAL-ranked bin panels A-F from a single h5ad dataset."
    )
    parser.add_argument("--config", required=True, help="YAML config file.")
    parser.add_argument("--h5ad", default=None, help="Optional override for input_h5ad in the config.")
    parser.add_argument("--output-dir", default=None, help="Optional override for output_dir in the config.")
    parser.add_argument("--patient-col", default=None, help="Optional override for obs_columns.patient_id.")
    parser.add_argument("--timepoint-col", default=None, help="Optional override for obs_columns.timepoint.")
    parser.add_argument("--response-col", default=None, help="Optional override for obs_columns.response.")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    input_h5ad = resolve_path(_coalesce(args.h5ad, config.get("input_h5ad")), base_dir=base_dir)
    output_dir = ensure_output_dir(resolve_path(_coalesce(args.output_dir, config.get("output_dir"), "figures/val_ranked_bin_panels"), base_dir=base_dir))

    ad = _load_anndata()
    if not input_h5ad.exists():
        raise FileNotFoundError(
            f"Input h5ad does not exist: {input_h5ad}\n"
            "For the example workflow, first run: python scripts/create_example_h5ad.py"
        )
    adata = ad.read_h5ad(input_h5ad)
    print(f"Loaded h5ad: {input_h5ad}")
    print(f"AnnData shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    cell_types = config.get("cell_types", list(CANONICAL_CELL_TYPES))
    obs_columns = config.get("obs_columns", {})
    sample_columns = SampleColumnConfig(
        patient_id=_coalesce(args.patient_col, obs_columns.get("patient_id")),
        timepoint=_coalesce(args.timepoint_col, obs_columns.get("timepoint")),
        response=_coalesce(args.response_col, obs_columns.get("response")),
        sample_id=obs_columns.get("sample_id"),
    )
    if not sample_columns.patient_id or not sample_columns.response:
        raise ValueError("Please provide obs_columns.patient_id and obs_columns.response in the config or CLI.")

    val_config = config.get("val", {})
    val_mode = val_config.get("mode", "marker_ensemble")
    if val_mode == "marker_ensemble":
        print("Computing VAL table with built-in marker ensemble mode.")
        val_table = build_val_table_from_marker_ensemble(
            adata,
            marker_sets=val_config.get("marker_sets"),
            layer=val_config.get("layer"),
            use_raw=bool(val_config.get("use_raw", False)),
            expression_transform=val_config.get("expression_transform", "auto_log1p"),
            min_markers_per_cell_type=int(val_config.get("min_markers_per_cell_type", 2)),
        )
    elif val_mode == "obs_annotations":
        print("Computing VAL table from annotation columns already present in adata.obs.")
        val_table = build_val_table_from_obs(
            adata,
            annotation_columns=val_config.get("annotation_columns", []),
            confidence_columns=val_config.get("confidence_columns", {}),
            label_map=val_config.get("label_map", {}),
        )
    else:
        raise ValueError("val.mode must be 'marker_ensemble' or 'obs_annotations'.")

    val_table_path = output_dir / "val_cell_table.csv"
    val_table.to_csv(val_table_path, index=False)
    print(f"Saved VAL cell table: {val_table_path}")

    binning_config = config.get("binning", {})
    n_bins = int(binning_config.get("n_bins", 10))
    denominator = binning_config.get("denominator", "all_cells")
    binned = assign_ranked_bins(val_table, cell_types=cell_types, n_bins=n_bins)
    binned_path = output_dir / "ranked_cell_bins.csv"
    binned.to_csv(binned_path, index=False)
    print(f"Saved ranked cell bins: {binned_path}")

    abundance = compute_ranked_bin_abundance(
        adata,
        binned,
        sample_columns=sample_columns,
        cell_types=cell_types,
        n_bins=n_bins,
        denominator=denominator,
    )
    abundance_path = output_dir / "ranked_bin_abundance_by_sample.csv"
    abundance.to_csv(abundance_path, index=False)
    print(f"Saved ranked-bin abundance table: {abundance_path}")

    response_config = config.get("response_groups", {})
    positive_response = response_config.get("positive", "R")
    negative_response = response_config.get("negative", "NR")
    stats_config = config.get("statistics", {})
    association = compute_ranked_bin_associations(
        abundance,
        positive_response=positive_response,
        negative_response=negative_response,
        fdr_scope=stats_config.get("fdr_scope", "global"),
    )
    association_path = output_dir / "ranked_bin_association_results.csv"
    association.to_csv(association_path, index=False)
    print(f"Saved ranked-bin association results: {association_path}")

    plot_config = config.get("plot", {})
    plot_ranked_bin_panels(
        association,
        output_dir=output_dir,
        cell_types=cell_types,
        fdr_threshold=float(plot_config.get("fdr_threshold", 0.05)),
        dpi=int(plot_config.get("dpi", 300)),
        basename=plot_config.get("basename", "ranked_bin_panels_A_to_F"),
    )
    print(f"Saved panels A-F to: {output_dir}")

    summary = (
        association.sort_values(["consensus_cell_type", "q_value"], ascending=[True, True])
        .groupby("consensus_cell_type", as_index=False)
        .first()[["consensus_cell_type", "bin_label", "bin_end_percent", "q_value", "neg_log10_q_value"]]
    )
    print("\nBest-ranked bin per cell type:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

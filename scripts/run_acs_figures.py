#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from val_toolkit.acs import (
    add_acs_rank_columns,
    annotation_table_from_csv,
    annotation_table_from_obs,
    compute_acs_cell_table,
    parse_annotation_methods,
)
from val_toolkit.annotation_harmonization import CANONICAL_ACS_CELL_TYPES
from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.val_panels import (
    SampleColumnConfig,
    assign_ranked_bins,
    build_val_table_from_marker_ensemble,
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
        description=(
            "Generate ACS figures A-F from a single h5ad dataset. ACS ranking follows the "
            "manuscript definition: VAL descending, then mean confidence among methods that "
            "support the consensus label."
        )
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
    output_dir = ensure_output_dir(
        resolve_path(_coalesce(args.output_dir, config.get("output_dir"), "figures/acs_figures"), base_dir=base_dir)
    )

    ad = _load_anndata()
    if not input_h5ad.exists():
        raise FileNotFoundError(
            f"Input h5ad does not exist: {input_h5ad}\n"
            "For the example workflow, first run: python scripts/create_example_h5ad.py"
        )
    adata = ad.read_h5ad(input_h5ad)
    print(f"Loaded h5ad: {input_h5ad}")
    print(f"AnnData shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    cell_types = config.get("major_cell_types", config.get("cell_types", list(CANONICAL_ACS_CELL_TYPES)))
    obs_columns = config.get("obs_columns", {})
    sample_columns = SampleColumnConfig(
        patient_id=_coalesce(args.patient_col, obs_columns.get("patient_id")),
        timepoint=_coalesce(args.timepoint_col, obs_columns.get("timepoint")),
        response=_coalesce(args.response_col, obs_columns.get("response")),
        sample_id=obs_columns.get("sample_id"),
    )
    if not sample_columns.patient_id or not sample_columns.response:
        raise ValueError("Please provide obs_columns.patient_id and obs_columns.response in the config or CLI.")

    annotation_source = config.get("annotation_source", {})
    mode = annotation_source.get("mode", "obs_columns")
    acs_config = config.get("acs", {})

    if mode == "obs_columns":
        methods = parse_annotation_methods(config.get("annotation_methods", []))
        print("Computing manuscript-style ACS from annotation columns in adata.obs.")
        annotation_table = annotation_table_from_obs(adata, methods)
        acs_cell_table = compute_acs_cell_table(
            annotation_table,
            methods=methods,
            label_mapping=config.get("label_mapping", {}),
            major_cell_types=cell_types,
        )
    elif mode == "annotation_csv":
        methods = parse_annotation_methods(config.get("annotation_methods", []))
        csv_path = resolve_path(annotation_source.get("path"), base_dir=base_dir)
        cell_id_column = annotation_source.get("cell_id_column", "cell_id")
        print(f"Computing manuscript-style ACS from external annotation CSV: {csv_path}")
        annotation_table = annotation_table_from_csv(csv_path, methods=methods, cell_id_column=cell_id_column)
        acs_cell_table = compute_acs_cell_table(
            annotation_table,
            methods=methods,
            label_mapping=config.get("label_mapping", {}),
            major_cell_types=cell_types,
        )
    elif mode == "marker_ensemble":
        # This is a self-contained fallback. It is ACS-like, but not full manuscript-style ACS
        # because the votes come from marker-scoring strategies rather than external annotation methods.
        marker_config = acs_config.get("marker_ensemble", {})
        print("Computing ACS-like table with built-in marker ensemble mode.")
        acs_cell_table = build_val_table_from_marker_ensemble(
            adata,
            marker_sets=marker_config.get("marker_sets"),
            layer=marker_config.get("layer"),
            use_raw=bool(marker_config.get("use_raw", False)),
            expression_transform=marker_config.get("expression_transform", "auto_log1p"),
            min_markers_per_cell_type=int(marker_config.get("min_markers_per_cell_type", 2)),
        )
        acs_cell_table["VAL"] = acs_cell_table["VAL_votes"]
        acs_cell_table["ACS_rank_score"] = acs_cell_table["VAL_rank_score"]
    else:
        raise ValueError("annotation_source.mode must be 'obs_columns', 'annotation_csv', or 'marker_ensemble'.")

    acs_cell_table = add_acs_rank_columns(acs_cell_table, major_cell_types=cell_types)

    acs_table_path = output_dir / "acs_cell_table.csv"
    acs_cell_table.to_csv(acs_table_path, index=False)
    print(f"Saved ACS cell table: {acs_table_path}")

    binning_config = config.get("binning", {})
    n_bins = int(binning_config.get("n_bins", 10))
    denominator = binning_config.get("denominator", "all_cells")
    binned = assign_ranked_bins(acs_cell_table, cell_types=cell_types, n_bins=n_bins)
    binned = binned.rename(
        columns={
            "rank_within_cell_type": "acs_rank_within_cell_type_from_bins",
            "rank_percentile": "acs_percentile_from_bins",
        }
    )
    binned_path = output_dir / "acs_ranked_cell_bins.csv"
    binned.to_csv(binned_path, index=False)
    print(f"Saved ACS-ranked cell bins: {binned_path}")

    abundance = compute_ranked_bin_abundance(
        adata,
        binned,
        sample_columns=sample_columns,
        cell_types=cell_types,
        n_bins=n_bins,
        denominator=denominator,
    )
    abundance_path = output_dir / "acs_ranked_bin_abundance_by_sample.csv"
    abundance.to_csv(abundance_path, index=False)
    print(f"Saved ACS-ranked bin abundance table: {abundance_path}")

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
    association_path = output_dir / "acs_ranked_bin_association_results.csv"
    association.to_csv(association_path, index=False)
    print(f"Saved ACS-ranked bin association results: {association_path}")

    plot_config = config.get("plot", {})
    basename = plot_config.get("basename", "acs_figures_A_to_F")
    plot_ranked_bin_panels(
        association,
        output_dir=output_dir,
        cell_types=cell_types,
        fdr_threshold=float(plot_config.get("fdr_threshold", 0.05)),
        dpi=int(plot_config.get("dpi", 300)),
        basename=basename,
        x_axis_prefix="ACS-ranked",
    )
    print(f"Saved ACS figures A-F to: {output_dir}")

    if not association.empty:
        summary = (
            association.sort_values(["consensus_cell_type", "q_value"], ascending=[True, True])
            .groupby("consensus_cell_type", as_index=False)
            .first()[["consensus_cell_type", "bin_label", "bin_end_percent", "q_value", "neg_log10_q_value"]]
        )
        print("\nBest ACS-ranked bin per cell type:")
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from val_toolkit.config import read_yaml, resolve_path
from val_toolkit.validation import ensure_output_dir


def _coalesce(cli_value: Any, config_value: Any, default: Any = None) -> Any:
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return default


def _repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def _write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)


def _run_step(command: list[str], label: str) -> None:
    print("\n" + "=" * 80)
    print(label)
    print("=" * 80)
    print("Running:", " ".join(command))
    subprocess.run(command, check=True)


def _any_backend_enabled(config: dict[str, Any]) -> bool:
    backends = config.get("annotation_backends", {}) or {}
    return any(bool(section.get("enabled", False)) for section in backends.values() if isinstance(section, dict))


def _split_column_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, float) and pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(",") if x.strip() and x.strip().lower() != "none"]


def _backend_method_name(backend_name: str, label_column: str) -> str:
    base = label_column
    if base.endswith("_label"):
        base = base[: -len("_label")]
    if backend_name.lower() == "singler" and base.startswith("singler_"):
        return "SingleR " + base[len("singler_") :].replace("_", " ").upper()
    if backend_name.lower() == "celltypist":
        return "CellTypist"
    return base.replace("_", " ").title()


def _augment_annotation_methods_from_status(
    annotation_methods: list[dict[str, Any]],
    status_csv: Path,
) -> list[dict[str, Any]]:
    if not status_csv.exists():
        return annotation_methods
    existing_label_cols = {str(m.get("label_column")) for m in annotation_methods if m.get("label_column")}
    status = pd.read_csv(status_csv)
    out = list(annotation_methods)
    for _, row in status.iterrows():
        if str(row.get("status", "")).lower() != "completed":
            continue
        label_cols = _split_column_list(row.get("label_column"))
        confidence_cols = _split_column_list(row.get("confidence_column"))
        if len(label_cols) != len(confidence_cols):
            raise ValueError(
                f"Backend status row for {row.get('name')} has mismatched label/confidence columns: "
                f"{label_cols} vs {confidence_cols}"
            )
        for label_col, confidence_col in zip(label_cols, confidence_cols):
            if label_col in existing_label_cols:
                continue
            out.append(
                {
                    "name": _backend_method_name(str(row.get("name", "Backend")), label_col),
                    "label_column": label_col,
                    "confidence_column": confidence_col,
                }
            )
            existing_label_cols.add(label_col)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the one-command ACS pipeline: optional annotation backends followed by "
            "ACS-ranked bin figures A-F."
        )
    )
    parser.add_argument("--config", required=True, help="YAML config file for the full ACS pipeline.")
    parser.add_argument("--h5ad", default=None, help="Optional override for input_h5ad in the config.")
    parser.add_argument("--output-dir", default=None, help="Optional override for output_dir in the config.")
    parser.add_argument(
        "--skip-backends",
        action="store_true",
        help="Skip annotation backend execution and run ACS directly on the input h5ad/annotation columns.",
    )
    args = parser.parse_args()

    repo_root = _repo_root_from_script()
    config_path = Path(args.config).resolve()
    config = read_yaml(config_path)
    base_dir = config_path.parent.parent

    input_h5ad = resolve_path(_coalesce(args.h5ad, config.get("input_h5ad")), base_dir=base_dir)
    output_dir = ensure_output_dir(
        resolve_path(_coalesce(args.output_dir, config.get("output_dir"), "figures/acs_pipeline"), base_dir=base_dir)
    )

    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad does not exist: {input_h5ad}")

    print(f"ACS pipeline input h5ad: {input_h5ad}")
    print(f"ACS pipeline output directory: {output_dir}")

    pipeline_cfg = config.get("pipeline", {}) or {}
    run_backends = bool(pipeline_cfg.get("run_annotation_backends", True)) and not args.skip_backends
    backend_enabled = _any_backend_enabled(config)

    annotation_output_dir = ensure_output_dir(output_dir / str(pipeline_cfg.get("annotation_output_subdir", "annotation_backends")))
    acs_output_dir = ensure_output_dir(output_dir / str(pipeline_cfg.get("acs_output_subdir", "acs_figures")))
    generated_config_dir = ensure_output_dir(output_dir / str(pipeline_cfg.get("generated_config_subdir", "generated_configs")))

    annotated_h5ad = input_h5ad
    backend_outputs = config.get("backend_outputs", {}) or {}
    status_csv_name = backend_outputs.get("status_csv", "annotation_backend_status.csv")

    if run_backends and backend_enabled:
        annotated_h5ad_name = backend_outputs.get("annotated_h5ad", "annotation_backends_annotated.h5ad")

        backend_config = {
            "input_h5ad": str(input_h5ad),
            "output_dir": str(annotation_output_dir),
            "expression": config.get("expression", {}) or {},
            "annotation_backends": config.get("annotation_backends", {}) or {},
            "outputs": {
                "combined_annotation_csv": backend_outputs.get("combined_annotation_csv", "annotation_backend_columns.csv"),
                "status_csv": status_csv_name,
                "write_h5ad": bool(backend_outputs.get("write_h5ad", True)),
                "annotated_h5ad": annotated_h5ad_name,
            },
        }
        backend_config_path = generated_config_dir / "annotation_backends.generated.yaml"
        _write_yaml(backend_config_path, backend_config)

        _run_step(
            [sys.executable, str(repo_root / "scripts" / "run_annotation_backends.py"), "--config", str(backend_config_path)],
            "Step 1/2: annotation backends",
        )
        if backend_config["outputs"].get("write_h5ad", True):
            annotated_h5ad = annotation_output_dir / annotated_h5ad_name
            if not annotated_h5ad.exists():
                raise FileNotFoundError(
                    "Annotation backend step was expected to write an annotated h5ad, but it was not found: "
                    f"{annotated_h5ad}"
                )
    else:
        reason = "--skip-backends was supplied" if args.skip_backends else "no annotation backends are enabled"
        print(f"\nSkipping annotation backend step because {reason}.")
        print("ACS will run directly on annotation columns already present in the input h5ad or external CSV.")

    annotation_methods = list(config.get("annotation_methods", []))
    if bool(pipeline_cfg.get("auto_include_backend_methods", True)) and run_backends and backend_enabled:
        annotation_methods = _augment_annotation_methods_from_status(
            annotation_methods,
            annotation_output_dir / status_csv_name,
        )
        print("\nACS annotation methods selected:")
        for method in annotation_methods:
            print(f"  - {method.get('name')}: {method.get('label_column')} / {method.get('confidence_column')}")

    acs_config = {
        "input_h5ad": str(annotated_h5ad),
        "output_dir": str(acs_output_dir),
        "obs_columns": config.get("obs_columns", {}) or {},
        "response_groups": config.get("response_groups", {}) or {},
        "major_cell_types": config.get("major_cell_types", ["B", "CD4 T", "CD8 T", "DC", "Monocytes", "NK"]),
        "annotation_source": config.get("annotation_source", {"mode": "obs_columns"}) or {"mode": "obs_columns"},
        "annotation_methods": annotation_methods,
        "label_mapping": config.get("label_mapping", {}) or {},
        "acs": config.get("acs", {}) or {},
        "binning": config.get("binning", {}) or {},
        "statistics": config.get("statistics", {}) or {},
        "plot": config.get("plot", {}) or {},
    }
    acs_config_path = generated_config_dir / "acs_figures.generated.yaml"
    _write_yaml(acs_config_path, acs_config)

    _run_step(
        [sys.executable, str(repo_root / "scripts" / "run_acs_figures.py"), "--config", str(acs_config_path)],
        "Step 2/2: ACS figures A-F",
    )

    print("\nACS pipeline complete.")
    print(f"Annotation backend outputs: {annotation_output_dir}")
    print(f"ACS figure outputs: {acs_output_dir}")
    print(f"Generated config files: {generated_config_dir}")


if __name__ == "__main__":
    main()

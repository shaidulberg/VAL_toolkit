from pathlib import Path

import yaml


def test_acs_pipeline_config_has_one_command_sections():
    cfg = yaml.safe_load(Path("configs/acs_pipeline.example.yaml").read_text())
    assert cfg["input_h5ad"].endswith(".h5ad")
    assert cfg["pipeline"]["run_annotation_backends"] is True
    assert cfg["annotation_backends"]["celltypist"]["enabled"] is True
    assert cfg["annotation_methods"][0]["label_column"] == "celltypist_label"
    assert cfg["annotation_methods"][0]["confidence_column"] == "celltypist_confidence"


def test_acs_pipeline_script_exists_and_is_executable():
    path = Path("scripts/run_acs_pipeline.py")
    assert path.exists()
    assert "run_annotation_backends.py" in path.read_text()
    assert "run_acs_figures.py" in path.read_text()


import pandas as pd

from scripts.run_acs_pipeline import _augment_annotation_methods_from_status


def test_pipeline_auto_includes_completed_backend_methods(tmp_path):
    status_csv = tmp_path / "annotation_backend_status.csv"
    pd.DataFrame(
        [
            {
                "name": "SingleR",
                "enabled": True,
                "status": "completed",
                "message": "ok",
                "output_csv": "singler.csv",
                "label_column": "singler_hpca_label,singler_monaco_label",
                "confidence_column": "singler_hpca_confidence,singler_monaco_confidence",
            }
        ]
    ).to_csv(status_csv, index=False)

    methods = _augment_annotation_methods_from_status([], status_csv)

    assert [m["label_column"] for m in methods] == ["singler_hpca_label", "singler_monaco_label"]
    assert [m["confidence_column"] for m in methods] == [
        "singler_hpca_confidence",
        "singler_monaco_confidence",
    ]

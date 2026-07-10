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

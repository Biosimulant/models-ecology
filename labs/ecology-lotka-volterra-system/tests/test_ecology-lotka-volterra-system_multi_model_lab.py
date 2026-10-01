from __future__ import annotations

from pathlib import Path

import yaml


def test_lab_manifest_uses_embedded_visualisation_model():
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())
    model_aliases = [entry["alias"] for entry in manifest["models"]]
    assert "visualisation" in model_aliases
    assert all(str(entry["path"]).startswith("models/") for entry in manifest["models"])

    assert not any(str(entry["maps_to"]).startswith("visualisation.") for entry in manifest["io"]["outputs"])

    wiring_targets = {target for entry in manifest["wiring"] for target in entry["to"]}
    assert any(target.startswith("visualisation.") for target in wiring_targets)


def test_lab_package_preserves_component_license_notices():
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())
    licenses = {entry["identifier"]: entry for entry in manifest["licenses"]}
    assert licenses["MIT"]["scope"] == "models/core/src/lotka_volterra.py"
    for entry in licenses.values():
        assert (lab_dir / entry["notice"]).is_file()
    core = (lab_dir / licenses["MIT"]["scope"]).read_text()
    assert "SPDX-License-Identifier: MIT" in core
    assert "does not relicense" in (lab_dir / "NOTICE.md").read_text()

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


def test_lab_license_scopes_preserve_core_and_native_source():
    import hashlib
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())
    licenses = {entry["identifier"]: entry for entry in manifest["licenses"]}
    assert "SPDX-License-Identifier: MIT" in (lab_dir / licenses["MIT"]["scope"]).read_text()
    source = lab_dir / licenses["CC0-1.0"]["scope"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == "27ccb5637322e70715e0ac40f50438fcc4a2b3167ae62c882a671c35473e55b5"
    for name in ("LICENSE-MIT.txt", "LICENSE-CODE.txt", "LICENSE-CONTENT.txt", "ATTRIBUTION.md"):
        assert (lab_dir / name).is_file()
    assert "not a license statement embedded" in (lab_dir / "NOTICE.md").read_text()

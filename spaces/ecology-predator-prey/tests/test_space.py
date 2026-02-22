from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest
import yaml


def _load_space():
    return yaml.safe_load((Path(__file__).resolve().parents[1] / "space.yaml").read_text(encoding="utf-8"))


def _repo_root_map(current_repo_root: Path) -> dict[str, Path]:
    return {"Biosimulant/models-ecology": current_repo_root.resolve()}


def _resolve_model_manifest(repo_map: dict[str, Path], model_ref: dict) -> Path:
    repo_full_name = str(model_ref.get("repo") or model_ref.get("repo_full_name") or "").strip()
    manifest_rel = str(model_ref.get("manifest_path") or "").strip()
    return (repo_map[repo_full_name] / manifest_rel).resolve()


def _clear_module_cache(module_name: str) -> None:
    root = module_name.split(".", 1)[0]
    to_delete = [k for k in sys.modules if k == root or k.startswith(f"{root}.")]
    for k in to_delete:
        sys.modules.pop(k, None)


def test_space_schema_and_paths():
    s = _load_space()
    assert s["schema_version"] == "2.0"
    assert s["models"]
    assert "runtime" in s and "wiring" in s
    current_repo_root = Path(__file__).resolve().parents[3]
    repo_map = _repo_root_map(current_repo_root)
    for m in s["models"]:
        assert _resolve_model_manifest(repo_map, m).exists()


def test_wiring_alias_references():
    s = _load_space()
    aliases = {m["alias"] for m in s["models"]}
    for w in s["wiring"]:
        src_alias = w["from"].split(".", 1)[0]
        assert src_alias in aliases
        for dst in w.get("to", []):
            dst_alias = dst.split(".", 1)[0]
            assert dst_alias in aliases


def test_space_smoke_runs_if_bsim_available():
    current_repo_root = Path(__file__).resolve().parents[3]
    monorepo_root = current_repo_root.parents[1]
    bsim_src = monorepo_root / "bsim" / "src"
    if bsim_src.exists():
        sys.path.insert(0, str(bsim_src))
    biosim = pytest.importorskip("biosim")
    s = _load_space()
    repo_map = _repo_root_map(current_repo_root)

    world = biosim.BioWorld()
    wb = biosim.WiringBuilder(world)
    for m in s["models"]:
        manifest_path = _resolve_model_manifest(repo_map, m)
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        meta = manifest.get("biosim") or {}
        ep = meta["entrypoint"]
        module_name, class_name = ep.split(":", 1)
        model_dir = manifest_path.parent
        sys.path.insert(0, str(model_dir))
        kwargs = dict(meta.get("init_kwargs") or {})
        kwargs.update(dict(m.get("parameters") or {}))
        _clear_module_cache(module_name)
        importlib.invalidate_caches()
        cls = getattr(importlib.import_module(module_name), class_name)
        wb.add(m["alias"], cls(**kwargs))

    for w in s["wiring"]:
        wb.connect(w["from"], w.get("to", []))
    wb.apply()
    tick_dt = float(s["runtime"]["tick_dt"])
    duration = min(float(s["runtime"]["duration"]), tick_dt * 20)
    world.run(duration=duration, tick_dt=tick_dt)
    visuals = world.collect_visuals()
    assert isinstance(visuals, list)

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import yaml


def _find_bsim_src(start: Path) -> Path | None:
    for parent in [start, *start.parents]:
        candidate = parent / "bsim-active" / "biosim" / "src"
        if (candidate / "biosim").is_dir():
            return candidate
    return None


def _ensure_paths() -> None:
    pack_root = Path(__file__).resolve().parents[1]
    if str(pack_root) not in sys.path:
        sys.path.insert(0, str(pack_root))

    bsim_src = _find_bsim_src(pack_root)
    if bsim_src is not None and str(bsim_src) not in sys.path:
        sys.path.insert(0, str(bsim_src))


def _load_module_class():
    _ensure_paths()
    manifest = Path(__file__).resolve().parents[1] / "model.yaml"
    data = yaml.safe_load(manifest.read_text(encoding="utf-8"))
    module_name, class_name = data["biosim"]["entrypoint"].split(":", 1)
    module = importlib.import_module(module_name)
    return getattr(module, class_name)


def _make_instance_and_advance():
    cls = _load_module_class()
    module = cls()
    t = float(getattr(module, "min_dt", 1.0) or 1.0)
    module.advance_to(t)
    return module, module.get_outputs()


def test_instantiation():
    cls = _load_module_class()
    module = cls()
    assert getattr(module, "min_dt", 0) > 0
    assert isinstance(module.inputs(), set)
    assert isinstance(module.outputs(), set)
    assert len(module.outputs()) > 0


def test_advance_produces_outputs():
    module, outputs = _make_instance_and_advance()
    assert isinstance(outputs, dict)
    assert set(outputs.keys()) == set(module.outputs())


def test_visualize_after_advance():
    module, _outputs = _make_instance_and_advance()
    visuals = module.visualize()
    assert isinstance(visuals, list)
    assert len(visuals) == 4

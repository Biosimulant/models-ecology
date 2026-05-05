from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import yaml


_MODEL_DIR = Path(__file__).resolve().parents[1]


def _find_bsim_src(start: Path) -> Path | None:
    for parent in [start, *start.parents]:
        for candidate in (parent / "biosim" / "src", parent / "bsim-active" / "biosim" / "src"):
            if (candidate / "biosim").is_dir():
                return candidate
    return None


_BSIM_SRC = _find_bsim_src(_MODEL_DIR)
if _BSIM_SRC is not None and str(_BSIM_SRC) not in sys.path:
    sys.path.insert(0, str(_BSIM_SRC))

from biosim.signals import RecordSignal


def _load_model():
    manifest = yaml.safe_load((_MODEL_DIR / "model.yaml").read_text())
    entrypoint = manifest["biosim"]["entrypoint"]
    module_name, class_name = entrypoint.split(":")
    module_rel = Path(*module_name.split(".")).with_suffix(".py")
    wrapper_path = _MODEL_DIR / module_rel
    unique_name = f"done_visualisation__{_MODEL_DIR.parent.parent.name}"
    spec = importlib.util.spec_from_file_location(unique_name, wrapper_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[unique_name] = module
    spec.loader.exec_module(module)
    module_cls = getattr(module, class_name)
    return module_cls(**dict(manifest["biosim"].get("init_kwargs", {}))), float(manifest["biosim"]["communication_step"])


def _sample_payload(mode: str) -> dict:
    if mode == "lotka_volterra":
        return {"parameters": {"alpha": 1.1, "beta": 0.4, "gamma": 0.4, "delta": 0.1, "prey_name": "Prey", "predator_name": "Predator"}, "prey_extinction_time": None, "predator_extinction_time": None, "history": [{"t": 0.0, "prey": 10.0, "predator": 5.0, "invariant": 1.0, "drift": 0.0}, {"t": 1.0, "prey": 8.0, "predator": 5.8, "invariant": 1.0, "drift": 0.0}]}
    if mode == "pfeiffer":
        return {"history": [{"t": 0.0, "S": 5.0, "N1": 2.0, "N2": 1.0}, {"t": 1.0, "S": 4.0, "N1": 3.0, "N2": 1.5}]}
    if mode == "turner":
        return {"history": [{"t": 0.0, "Population_of_Eggs": 10.0, "Population_of_Larvae": 4.0, "Population_of_Pupae": 1.0}, {"t": 1.0, "Population_of_Eggs": 9.0, "Population_of_Larvae": 5.0, "Population_of_Pupae": 2.0}]}
    if mode == "leibovich":
        return {"species_count": 3, "history": [{"t": 0.0, "species_1": 10.0, "species_2": 8.0, "species_3": 5.0, "total_abundance": 23.0, "richness": 3.0, "shannon_diversity": 1.05, "dominant_species_index": 1.0}, {"t": 1.0, "species_1": 12.0, "species_2": 7.0, "species_3": 4.0, "total_abundance": 23.0, "richness": 3.0, "shannon_diversity": 1.0, "dominant_species_index": 1.0}]}
    if mode == "gene_drive":
        return {"initial_total": 100.0, "history": [{"t": 1.0, "total_adults": 90.0, "adult_females": 40.0, "adult_males": 50.0, "drive_frequency": 0.4, "resistance_frequency": 0.1, "male_fraction": 0.56, "suppression_ratio": 0.1}, {"t": 2.0, "total_adults": 80.0, "adult_females": 33.0, "adult_males": 47.0, "drive_frequency": 0.55, "resistance_frequency": 0.15, "male_fraction": 0.5875, "suppression_ratio": 0.2}]}
    return {"time_unit": "day", "prey_label": "Prey", "predator_label": "Predator", "equilibrium_summary": {"prey_equilibrium": 12.0, "predator_equilibrium": 4.0}, "stability_summary": {"regime": "stable"}, "extinction_risk": {"prey": 0.1, "predator": 0.2}, "history": [{"t": 0.0, "prey": 20.0, "predator": 4.0}, {"t": 1.0, "prey": 18.0, "predator": 4.5}]}


def test_visualisation_model_stays_internal_and_renders_visuals():
    module, step = _load_model()
    alias = module.source_alias
    spec = module.inputs()[f"{alias}_visualisation_payload"]
    payload = _sample_payload(module.mode)
    module.set_inputs({f"{alias}_visualisation_payload": RecordSignal(source="test", name=f"{alias}_visualisation_payload", value={"payload": payload}, emitted_at=step, spec=spec)})
    module.advance_window(0.0, step)

    assert module.outputs() == {}
    assert module.get_outputs() == {}
    visuals = module.visualize()
    assert isinstance(visuals, list) and visuals
    renders = {visual["render"] for visual in visuals}
    assert "timeseries" in renders
    assert "table" in renders

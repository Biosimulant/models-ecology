from __future__ import annotations
from biosim.signals import (AcceptedSignalProfile, ArraySignal, BioSignal, EventSignal, RecordSignal, ScalarSignal, SignalSpec)

import sys
from pathlib import Path


MODEL_ROOT = Path(__file__).resolve().parents[1]
MONOREPO_ROOT = MODEL_ROOT.parents[3]
BSIM_SRC = MONOREPO_ROOT / "bsim-active" / "biosim" / "src"

for path in (str(MODEL_ROOT), str(BSIM_SRC)):
    if path not in sys.path:
        sys.path.insert(0, path)


from src.turner2015_human_mosquito_elp_model_biomd0000000922_model import (  # noqa: E402
    Turner2015HumanMosquitoElpModelBiomd0000000922Model,
)


class FakeRunner:
    def __init__(self) -> None:
        self._species = ["Population_of_Eggs", "Population_of_Larvae", "Population_of_Pupae"]

    def getFloatingSpeciesIds(self):
        return list(self._species)

    def simulate(self, start: float, end: float, steps: int, selections=None):
        points = []
        for idx in range(steps):
            t = start + (end - start) * idx / (steps - 1)
            points.append([t, 1000.0 - 5.0 * t, 700.0 + 2.0 * t, 200.0 + 3.0 * t])
        return points

    def reset(self) -> None:
        return None


def test_history_and_outputs_accumulate(monkeypatch) -> None:
    monkeypatch.setattr(
        Turner2015HumanMosquitoElpModelBiomd0000000922Model,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015HumanMosquitoElpModelBiomd0000000922Model(integration_step=0.5)
    module.advance_window(0.0, 2.0)
    outputs = module.get_outputs()

    assert set(outputs) == {"life_stage_state", "population_metrics"}
    assert len(module._history) > 1
    assert outputs["life_stage_state"].value["total_immature_population"] > 0.0


def test_visuals_are_multi_point(monkeypatch) -> None:
    monkeypatch.setattr(
        Turner2015HumanMosquitoElpModelBiomd0000000922Model,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015HumanMosquitoElpModelBiomd0000000922Model(integration_step=0.5)
    module.advance_window(0.0, 3.0)
    visuals = module.visualize()

    assert isinstance(visuals, list)
    assert [visual["render"] for visual in visuals] == ["timeseries", "timeseries", "table"]
    for visual in visuals[:2]:
        for series in visual["data"]["series"]:
            assert len(series["points"]) > 1


def test_stage_fractions_sum_to_one(monkeypatch) -> None:
    monkeypatch.setattr(
        Turner2015HumanMosquitoElpModelBiomd0000000922Model,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015HumanMosquitoElpModelBiomd0000000922Model(integration_step=0.5)
    module.advance_window(0.0, 1.0)
    metrics = module.get_outputs()["population_metrics"].value

    total = metrics["egg_fraction"] + metrics["larval_fraction"] + metrics["pupal_fraction"]
    assert abs(total - 1.0) < 1e-9

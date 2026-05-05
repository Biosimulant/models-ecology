from __future__ import annotations

from src.turner2015_mosquito_life_stages import Turner2015MosquitoLifeStagesModel


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
        Turner2015MosquitoLifeStagesModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015MosquitoLifeStagesModel(integration_step=0.5)
    module.advance_window(0.0, 2.0)
    outputs = module.get_outputs()

    assert set(outputs) == {"life_stage_state", "population_metrics", "visualisation_payload"}
    assert len(module._history) > 1
    assert outputs["life_stage_state"].value["total_immature_population"] > 0.0


def test_visualisation_payload_is_multi_point(monkeypatch) -> None:
    monkeypatch.setattr(
        Turner2015MosquitoLifeStagesModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015MosquitoLifeStagesModel(integration_step=0.5)
    module.advance_window(0.0, 3.0)
    payload = module.get_outputs()["visualisation_payload"].value["payload"]
    assert len(payload["history"]) > 1


def test_stage_fractions_sum_to_one(monkeypatch) -> None:
    monkeypatch.setattr(
        Turner2015MosquitoLifeStagesModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Turner2015MosquitoLifeStagesModel(integration_step=0.5)
    module.advance_window(0.0, 1.0)
    metrics = module.get_outputs()["population_metrics"].value

    total = metrics["egg_fraction"] + metrics["larval_fraction"] + metrics["pupal_fraction"]
    assert abs(total - 1.0) < 1e-9

from __future__ import annotations

from src.pfeiffer2001_atp_cooperation import Pfeiffer2001AtpCooperationModel


class FakeRunner:
    def __init__(self) -> None:
        self._species = ["S", "N1", "N2"]

    def getFloatingSpeciesIds(self):
        return list(self._species)

    def simulate(self, start: float, end: float, steps: int, selections=None):
        points = []
        for idx in range(steps):
            t = start + (end - start) * idx / (steps - 1)
            points.append([t, max(0.0, 1.0 - 0.1 * t), 100.0 + 2.0 * t, 10.0 + 5.0 * t])
        return points

    def reset(self) -> None:
        return None


def test_history_accumulates_and_outputs_are_domain_specific(monkeypatch) -> None:
    monkeypatch.setattr(
        Pfeiffer2001AtpCooperationModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Pfeiffer2001AtpCooperationModel(integration_step=0.5)
    module.advance_window(0.0, 2.0)
    outputs = module.get_outputs()

    assert set(outputs) == {"resource_state", "community_state", "cooperation_metrics"}
    assert len(module._history) > 1
    assert outputs["community_state"].value["total_population"] > 0.0


def test_visuals_have_real_trajectories(monkeypatch) -> None:
    monkeypatch.setattr(
        Pfeiffer2001AtpCooperationModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Pfeiffer2001AtpCooperationModel(integration_step=0.5)
    module.advance_window(0.0, 3.0)
    visuals = module.visualize()

    assert isinstance(visuals, list)
    assert [visual["render"] for visual in visuals] == ["timeseries", "timeseries", "table"]
    for visual in visuals[:2]:
        for series in visual["data"]["series"]:
            assert len(series["points"]) > 1


def test_fraction_metrics_sum_to_one(monkeypatch) -> None:
    monkeypatch.setattr(
        Pfeiffer2001AtpCooperationModel,
        "_build_runner",
        lambda self: FakeRunner(),
    )

    module = Pfeiffer2001AtpCooperationModel(integration_step=0.5)
    module.advance_window(0.0, 1.0)
    metrics = module.get_outputs()["cooperation_metrics"].value

    assert abs(metrics["high_yield_fraction"] + metrics["low_yield_fraction"] - 1.0) < 1e-9

from __future__ import annotations
from biosim.signals import (AcceptedSignalProfile, ArraySignal, BioSignal, EventSignal, RecordSignal, ScalarSignal, SignalSpec)

import sys
from pathlib import Path

import pytest


MODEL_ROOT = Path(__file__).resolve().parents[1]
MONOREPO_ROOT = MODEL_ROOT.parents[3]
BSIM_SRC = MONOREPO_ROOT / "bsim-active" / "biosim" / "src"

for path in (str(MODEL_ROOT), str(BSIM_SRC)):
    if path not in sys.path:
        sys.path.insert(0, path)


from src.geci2022_model2301120001_model import Geci2022Model2301120001Model  # noqa: E402


def test_outputs_and_history_accumulate() -> None:
    module = Geci2022Model2301120001Model(integration_step=0.2)
    module.advance_window(0.0, 2.0)

    outputs = module.get_outputs()
    assert set(outputs) == {"population_state", "gene_drive_metrics"}
    assert len(module._history) == 10
    assert outputs["population_state"].value["adult_females"] < outputs["population_state"].value["total_adults"]


def test_drive_metrics_stay_in_physical_bounds() -> None:
    module = Geci2022Model2301120001Model(integration_step=0.1)
    module.advance_window(0.0, 5.0)

    metrics = module.get_outputs()["gene_drive_metrics"].value
    assert 0.0 <= metrics["drive_frequency"] <= 1.0
    assert 0.0 <= metrics["resistance_frequency"] <= 1.0
    assert 0.5 <= metrics["male_fraction"] <= 0.98


def test_visuals_use_multi_point_trajectories() -> None:
    module = Geci2022Model2301120001Model(integration_step=0.25)
    module.advance_window(0.0, 5.0)
    visuals = module.visualize()

    assert isinstance(visuals, list)
    assert [visual["render"] for visual in visuals] == ["timeseries", "timeseries", "table"]
    for series in visuals[0]["data"]["series"] + visuals[1]["data"]["series"]:
        assert len(series["points"]) > 1


def test_suppression_increases_under_stronger_drive() -> None:
    weak = Geci2022Model2301120001Model(suppression_strength=0.2, integration_step=0.1)
    strong = Geci2022Model2301120001Model(suppression_strength=1.2, integration_step=0.1)

    weak.advance_window(0.0, 6.0)
    strong.advance_window(0.0, 6.0)

    weak_population = weak.get_outputs()["population_state"].value["total_adults"]
    strong_population = strong.get_outputs()["population_state"].value["total_adults"]
    assert strong_population < weak_population

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml


def _reference_rk4(alpha: float, beta: float, gamma: float, delta: float, prey: float, predator: float, dt: float, steps: int) -> tuple[float, float]:
    def rhs(x: float, y: float) -> tuple[float, float]:
        return alpha * x - beta * x * y, delta * x * y - gamma * y

    for _ in range(steps):
        k1x, k1y = rhs(prey, predator)
        k2x, k2y = rhs(prey + 0.5 * dt * k1x, predator + 0.5 * dt * k1y)
        k3x, k3y = rhs(prey + 0.5 * dt * k2x, predator + 0.5 * dt * k2y)
        k4x, k4y = rhs(prey + dt * k3x, predator + dt * k3y)
        prey += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        predator += (dt / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
    return prey, predator


def test_instantiation(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem()
    assert module.integration_step > 0
    assert set(module.inputs()) == {
        "prey_initial_population",
        "predator_initial_population",
        "prey_growth_rate",
        "predation_rate",
        "predator_mortality_rate",
        "predator_reproduction_rate",
    }
    outputs = module.outputs()
    assert set(outputs) == {"prey_population_state", "predator_population_state", "visualisation_payload"}
    assert outputs["prey_population_state"].emitted_unit == "count"
    assert outputs["predator_population_state"].emitted_unit == "count"


def test_prey_grows_without_predators(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(predator_initial=0.0, prey_initial=10.0, integration_step=0.1)
    module.advance_window(0.0, 2.0)
    outputs = module.get_outputs()
    assert outputs["prey_population_state"].value["count"] > 10.0
    assert outputs["predator_population_state"].value["count"] == 0.0


def test_predators_decay_without_prey(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(prey_initial=0.0, predator_initial=5.0, integration_step=0.1)
    module.advance_window(0.0, 2.0)
    outputs = module.get_outputs()
    assert outputs["prey_population_state"].value["count"] == 0.0
    assert outputs["predator_population_state"].value["count"] < 5.0


def test_reference_trajectory_matches_small_step_reference(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    alpha, beta, gamma, delta = 1.1, 0.4, 0.4, 0.1
    prey_initial, predator_initial = 10.0, 5.0
    module = LotkaVolterraSystem(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        delta=delta,
        prey_initial=prey_initial,
        predator_initial=predator_initial,
        integration_step=0.1,
    )
    module.advance_window(0.0, 1.0)
    outputs = module.get_outputs()
    prey = outputs["prey_population_state"].value["count"]
    predator = outputs["predator_population_state"].value["count"]
    ref_prey, ref_predator = _reference_rk4(alpha, beta, gamma, delta, prey_initial, predator_initial, 0.001, 1000)
    assert prey == pytest.approx(ref_prey, rel=0.0, abs=2e-4)
    assert predator == pytest.approx(ref_predator, rel=0.0, abs=2e-4)


def test_visualisation_payload_contains_expected_lv_fields(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(integration_step=0.05)
    module.advance_window(0.0, 20.0)
    payload = module.get_outputs()["visualisation_payload"].value["payload"]

    assert set(payload) == {"parameters", "prey_extinction_time", "predator_extinction_time", "point"}
    assert payload["point"]["t"] == pytest.approx(20.0)
    assert payload["parameters"]["alpha"] == 1.1
    assert payload["parameters"]["prey_name"] == "Prey"


def test_audit_drift_stays_bounded(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(integration_step=0.05)
    module.advance_window(0.0, 20.0)
    payload = module._visualisation_payload()
    max_abs_drift = max(abs(point["drift"]) for point in payload["history"])
    assert max_abs_drift < 1e-5


def test_populations_remain_finite_and_nonnegative(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(integration_step=0.05)
    module.advance_window(0.0, 30.0)
    outputs = module.get_outputs()
    prey_state = outputs["prey_population_state"].value
    predator_state = outputs["predator_population_state"].value
    assert prey_state["role"] == "prey"
    assert prey_state["label"] == "Prey"
    assert predator_state["role"] == "predator"
    assert predator_state["label"] == "Predator"
    assert outputs["prey_population_state"].spec.emitted_unit == "count"
    assert outputs["predator_population_state"].spec.emitted_unit == "count"
    assert math.isfinite(prey_state["count"]) and prey_state["count"] >= 0.0
    assert math.isfinite(predator_state["count"]) and predator_state["count"] >= 0.0


def test_lotka_volterra_equations_remain_unchanged(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(alpha=1.1, beta=0.4, gamma=0.4, delta=0.1)
    dprey, dpredator = module._rhs(prey=10.0, predator=5.0)
    assert dprey == pytest.approx(1.1 * 10.0 - 0.4 * 10.0 * 5.0)
    assert dpredator == pytest.approx(0.1 * 10.0 * 5.0 - 0.4 * 5.0)


def test_lab_level_io_maps_to_internal_model_outputs():
    lab_manifest = Path(__file__).resolve().parents[3] / "lab.yaml"
    data = yaml.safe_load(lab_manifest.read_text(encoding="utf-8"))
    output_maps = {entry["name"]: entry["maps_to"] for entry in data["io"]["outputs"]}
    assert output_maps["prey_population_state"] == "ecology_lotka_volterra_system.prey_population_state"
    assert output_maps["predator_population_state"] == "ecology_lotka_volterra_system.predator_population_state"


def test_rabbit_fox_labels_do_not_leak_into_wire_ports(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem()
    port_names = set(module.inputs()) | set(module.outputs())
    assert all("rabbit" not in name.lower() for name in port_names)
    assert all("fox" not in name.lower() for name in port_names)

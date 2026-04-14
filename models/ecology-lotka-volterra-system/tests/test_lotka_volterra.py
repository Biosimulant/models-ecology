from __future__ import annotations

import math

import pytest


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
    assert module.min_dt > 0
    assert module.inputs() == set()
    assert module.outputs() == {"prey_population_state", "predator_population_state"}


def test_prey_grows_without_predators(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(predator_initial=0.0, prey_initial=10.0, min_dt=0.1)
    module.advance_to(2.0)
    outputs = module.get_outputs()
    assert outputs["prey_population_state"].value["count"] > 10.0
    assert outputs["predator_population_state"].value["count"] == 0.0


def test_predators_decay_without_prey(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(prey_initial=0.0, predator_initial=5.0, min_dt=0.1)
    module.advance_to(2.0)
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
        min_dt=0.1,
    )
    module.advance_to(1.0)
    outputs = module.get_outputs()
    prey = outputs["prey_population_state"].value["count"]
    predator = outputs["predator_population_state"].value["count"]
    ref_prey, ref_predator = _reference_rk4(alpha, beta, gamma, delta, prey_initial, predator_initial, 0.001, 1000)
    assert prey == pytest.approx(ref_prey, rel=0.0, abs=2e-4)
    assert predator == pytest.approx(ref_predator, rel=0.0, abs=2e-4)


def test_visualize_returns_expected_lv_visuals(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(min_dt=0.05)
    module.advance_to(20.0)
    visuals = module.visualize()

    assert isinstance(visuals, list)
    assert len(visuals) == 4
    assert [visual["render"] for visual in visuals] == ["timeseries", "image", "table", "timeseries"]
    assert "equilibrium" in module._build_phase_svg()


def test_audit_drift_stays_bounded(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(min_dt=0.05)
    module.advance_to(20.0)
    visuals = module.visualize()
    audit_points = visuals[3]["data"]["series"][1]["points"]
    max_abs_drift = max(abs(point[1]) for point in audit_points)
    assert max_abs_drift < 1e-5


def test_populations_remain_finite_and_nonnegative(biosim):
    from src.lotka_volterra import LotkaVolterraSystem

    module = LotkaVolterraSystem(min_dt=0.05)
    module.advance_to(30.0)
    outputs = module.get_outputs()
    prey = outputs["prey_population_state"].value["count"]
    predator = outputs["predator_population_state"].value["count"]
    assert math.isfinite(prey) and prey >= 0.0
    assert math.isfinite(predator) and predator >= 0.0

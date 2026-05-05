from __future__ import annotations

import pytest

from src.rosenzweig_macarthur import RosenzweigMacArthurPredatorPreySystem


def run_days(module: RosenzweigMacArthurPredatorPreySystem, days: float) -> dict:
    step = 0.5
    target = module.get_state()["time"]
    for _ in range(int(days / step)):
        target += step
        module.advance_window(target - step, target)
    return module.get_outputs()


def value(outputs: dict, name: str) -> dict:
    return outputs[name].value


def test_carrying_capacity_bounds_prey_growth_without_predators():
    module = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=10,
        predator_initial_population=0,
        prey_carrying_capacity=80,
        prey_growth_rate=0.9,
        attack_rate=0.0,
    )
    outputs = run_days(module, 80)
    prey = value(outputs, "prey_population_state")["count"]
    assert prey <= 82
    assert prey > 70


def test_holling_type_ii_predation_saturates_at_high_prey_density():
    module = RosenzweigMacArthurPredatorPreySystem(attack_rate=0.08, handling_time=0.25)
    low = module._predation(1000, 10)
    high = module._predation(2000, 10)
    assert high > low
    assert high / low < 1.1


def test_zero_predator_population_stays_zero_unless_introduced():
    module = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=40,
        predator_initial_population=0,
        noise_strength=0,
        migration_rate=0,
    )
    outputs = run_days(module, 20)
    assert value(outputs, "predator_population_state")["count"] == pytest.approx(0.0)


def test_seasonality_changes_trajectories_reproducibly():
    baseline = RosenzweigMacArthurPredatorPreySystem(seasonal_amplitude=0.0)
    seasonal_a = RosenzweigMacArthurPredatorPreySystem(seasonal_amplitude=0.4, seasonal_period=30)
    seasonal_b = RosenzweigMacArthurPredatorPreySystem(seasonal_amplitude=0.4, seasonal_period=30)
    base_outputs = run_days(baseline, 17)
    seasonal_outputs_a = run_days(seasonal_a, 17)
    seasonal_outputs_b = run_days(seasonal_b, 17)
    assert value(seasonal_outputs_a, "prey_population_state")["count"] != pytest.approx(
        value(base_outputs, "prey_population_state")["count"]
    )
    assert value(seasonal_outputs_a, "prey_population_state")["count"] == pytest.approx(
        value(seasonal_outputs_b, "prey_population_state")["count"]
    )


def test_noise_is_reproducible_with_seed():
    first = RosenzweigMacArthurPredatorPreySystem(noise_strength=0.01, random_seed=42)
    second = RosenzweigMacArthurPredatorPreySystem(noise_strength=0.01, random_seed=42)
    outputs_a = run_days(first, 20)
    outputs_b = run_days(second, 20)
    assert value(outputs_a, "prey_population_state")["count"] == pytest.approx(
        value(outputs_b, "prey_population_state")["count"]
    )
    assert value(outputs_a, "predator_population_state")["count"] == pytest.approx(
        value(outputs_b, "predator_population_state")["count"]
    )


def test_food_limitation_lowers_effective_carrying_capacity():
    full = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=10,
        predator_initial_population=0,
        prey_carrying_capacity=100,
        food_resource_index=1.0,
        attack_rate=0.0,
    )
    limited = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=10,
        predator_initial_population=0,
        prey_carrying_capacity=100,
        food_resource_index=0.4,
        attack_rate=0.0,
    )
    full_outputs = run_days(full, 80)
    limited_outputs = run_days(limited, 80)
    assert value(limited_outputs, "prey_population_state")["count"] < value(
        full_outputs, "prey_population_state"
    )["count"]


def test_disease_reduces_susceptible_and_infected_compartments():
    diseased = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=100,
        predator_initial_population=0,
        prey_carrying_capacity=120,
        disease_enabled=True,
        prey_infected_initial=30,
        prey_transmission_rate=0.3,
        prey_recovery_rate=0.0,
        prey_disease_mortality_rate=0.25,
        attack_rate=0.0,
    )
    healthy = RosenzweigMacArthurPredatorPreySystem(
        prey_initial_population=100,
        predator_initial_population=0,
        prey_carrying_capacity=120,
        attack_rate=0.0,
    )
    diseased_outputs = run_days(diseased, 20)
    healthy_outputs = run_days(healthy, 20)
    assert sum(diseased._prey_infected) < 30
    assert value(diseased_outputs, "prey_population_state")["count"] < value(
        healthy_outputs, "prey_population_state"
    )["count"]


def test_migration_moves_population_between_patches_and_can_rescue_low_patch():
    module = RosenzweigMacArthurPredatorPreySystem(
        patch_count=2,
        patch_initial_prey=[100, 0],
        patch_initial_predators=[0, 0],
        predator_initial_population=0,
        attack_rate=0.0,
        migration_rate=0.2,
    )
    outputs = run_days(module, 5)
    patches = value(outputs, "prey_population_state")["patches"]
    assert patches[1] > 0
    assert patches[0] > patches[1]


def test_outputs_always_include_units_and_typed_state_records():
    module = RosenzweigMacArthurPredatorPreySystem()
    module.advance_window(0.0, 0.5)
    outputs = module.get_outputs()
    specs = module.outputs()
    assert specs["prey_population_state"].emitted_unit == "count"
    assert specs["predator_population_state"].emitted_unit == "count"
    assert outputs["prey_population_state"].spec.emitted_unit == "count"
    prey_state = value(outputs, "prey_population_state")
    assert prey_state["role"] == "prey"
    assert prey_state["label"] == "Prey"
    assert isinstance(prey_state["patches"], list)
    assert "equilibrium_summary" in outputs
    assert "scenario_summary" in outputs


def test_visualisation_payload_is_ecology_specific():
    module = RosenzweigMacArthurPredatorPreySystem(
        seasonal_amplitude=0.2,
        patch_count=2,
        patch_initial_prey=[80, 20],
        patch_initial_predators=[12, 4],
        migration_rate=0.05,
    )
    run_days(module, 10)
    payload = value(module.get_outputs(), "visualisation_payload")["payload"]
    assert payload["prey_label"] == "Prey"
    assert payload["predator_label"] == "Predator"
    assert len(payload["history"]) > 1


def test_visualisation_payload_spans_the_retained_run():
    module = RosenzweigMacArthurPredatorPreySystem()
    run_days(module, 200)

    prey_points = value(module.get_outputs(), "visualisation_payload")["payload"]["history"]

    assert prey_points[0]["t"] == pytest.approx(0.5)
    assert prey_points[-1]["t"] == pytest.approx(200.0)

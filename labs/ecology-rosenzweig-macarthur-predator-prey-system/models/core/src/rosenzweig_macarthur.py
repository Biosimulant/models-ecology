"""Rosenzweig-MacArthur predator-prey system for applied ecology labs."""
from __future__ import annotations

import base64
import math
import random
from typing import Any, Dict, List, Optional, Sequence, TYPE_CHECKING

from biosim import BioModule
from biosim.signals import AcceptedSignalProfile, BioSignal, RecordSignal, ScalarSignal, SignalSpec
from biosim.signals import coerce_float, unwrap_payload
from biosim.signals import make_signal as _make_signal

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


EPS = 1e-9
HISTORY_POINT_LIMIT = 10000
VISUAL_POINT_LIMIT = 900



class RosenzweigMacArthurPredatorPreySystem(BioModule):
    """Logistic prey growth with Holling type II predation and applied ecology mechanisms."""

    def __init__(
        self,
        prey_initial_population: float = 60.0,
        predator_initial_population: float = 12.0,
        prey_growth_rate: float = 0.8,
        prey_carrying_capacity: float = 120.0,
        attack_rate: float = 0.08,
        handling_time: float = 0.25,
        predator_mortality_rate: float = 0.25,
        predator_conversion_efficiency: float = 0.12,
        food_resource_index: float = 1.0,
        seasonal_amplitude: float = 0.0,
        seasonal_period: float = 365.0,
        seasonal_phase: float = 0.0,
        noise_strength: float = 0.0,
        random_seed: int | None = 1,
        age_structure_enabled: bool = False,
        juvenile_prey_fraction: float = 0.35,
        juvenile_predator_fraction: float = 0.25,
        prey_maturation_rate: float = 0.15,
        predator_maturation_rate: float = 0.08,
        juvenile_prey_mortality_rate: float = 0.02,
        juvenile_predator_mortality_rate: float = 0.03,
        disease_enabled: bool = False,
        prey_infected_initial: float = 0.0,
        predator_infected_initial: float = 0.0,
        prey_transmission_rate: float = 0.0,
        predator_transmission_rate: float = 0.0,
        prey_recovery_rate: float = 0.0,
        predator_recovery_rate: float = 0.0,
        prey_disease_mortality_rate: float = 0.0,
        predator_disease_mortality_rate: float = 0.0,
        patch_count: int = 1,
        migration_rate: float = 0.0,
        migration_matrix: Sequence[Sequence[float]] | None = None,
        patch_carrying_capacities: Sequence[float] | None = None,
        patch_food_resources: Sequence[float] | None = None,
        patch_initial_prey: Sequence[float] | None = None,
        patch_initial_predators: Sequence[float] | None = None,
        prey_extinction_threshold: float = 1.0,
        predator_extinction_threshold: float = 1.0,
        prey_label: str = "Prey",
        predator_label: str = "Predator",
        time_unit: str = "day",
        integration_step: float = 0.05,
        **_: Any,
    ) -> None:
        self.prey_growth_rate = max(0.0, float(prey_growth_rate))
        self.prey_carrying_capacity = max(EPS, float(prey_carrying_capacity))
        self.attack_rate = max(0.0, float(attack_rate))
        self.handling_time = max(0.0, float(handling_time))
        self.predator_mortality_rate = max(0.0, float(predator_mortality_rate))
        self.predator_conversion_efficiency = max(0.0, float(predator_conversion_efficiency))
        self.food_resource_index = max(0.0, float(food_resource_index))
        self.seasonal_amplitude = self._fraction(seasonal_amplitude, upper=0.95)
        self.seasonal_period = max(EPS, float(seasonal_period))
        self.seasonal_phase = float(seasonal_phase)
        self.noise_strength = max(0.0, float(noise_strength))
        self.random_seed = random_seed
        self._rng = random.Random(random_seed)

        self.age_structure_enabled = bool(age_structure_enabled)
        self.juvenile_prey_fraction = self._fraction(juvenile_prey_fraction)
        self.juvenile_predator_fraction = self._fraction(juvenile_predator_fraction)
        self.prey_maturation_rate = max(0.0, float(prey_maturation_rate))
        self.predator_maturation_rate = max(0.0, float(predator_maturation_rate))
        self.juvenile_prey_mortality_rate = max(0.0, float(juvenile_prey_mortality_rate))
        self.juvenile_predator_mortality_rate = max(0.0, float(juvenile_predator_mortality_rate))

        self.disease_enabled = bool(disease_enabled)
        self.prey_transmission_rate = max(0.0, float(prey_transmission_rate))
        self.predator_transmission_rate = max(0.0, float(predator_transmission_rate))
        self.prey_recovery_rate = max(0.0, float(prey_recovery_rate))
        self.predator_recovery_rate = max(0.0, float(predator_recovery_rate))
        self.prey_disease_mortality_rate = max(0.0, float(prey_disease_mortality_rate))
        self.predator_disease_mortality_rate = max(0.0, float(predator_disease_mortality_rate))

        self.patch_count = max(1, int(patch_count))
        self.migration_rate = max(0.0, float(migration_rate))
        self.migration_matrix = self._normalize_matrix(migration_matrix, self.patch_count)
        self.patch_carrying_capacities = self._expand_series(
            patch_carrying_capacities, self.patch_count, self.prey_carrying_capacity
        )
        self.patch_food_resources = self._expand_series(patch_food_resources, self.patch_count, 1.0)
        self.prey_initial = self._expand_series(patch_initial_prey, self.patch_count, prey_initial_population)
        self.predator_initial = self._expand_series(patch_initial_predators, self.patch_count, predator_initial_population)

        self.prey_extinction_threshold = max(0.0, float(prey_extinction_threshold))
        self.predator_extinction_threshold = max(0.0, float(predator_extinction_threshold))
        self.prey_label = prey_label or "Prey"
        self.predator_label = predator_label or "Predator"
        self.time_unit = time_unit or "day"
        self.integration_step = max(1e-4, float(integration_step))
        self.communication_step = getattr(self, "communication_step", self.integration_step)

        self._time = 0.0
        self._prey = [max(0.0, float(value)) for value in self.prey_initial]
        self._predator = [max(0.0, float(value)) for value in self.predator_initial]
        self._prey_infected = self._expand_series(None, self.patch_count, prey_infected_initial)
        self._predator_infected = self._expand_series(None, self.patch_count, predator_infected_initial)
        self._prey_susceptible = [max(0.0, total - infected) for total, infected in zip(self._prey, self._prey_infected)]
        self._predator_susceptible = [
            max(0.0, total - infected) for total, infected in zip(self._predator, self._predator_infected)
        ]
        self._prey_juvenile = [value * self.juvenile_prey_fraction for value in self._prey]
        self._prey_adult = [value - juvenile for value, juvenile in zip(self._prey, self._prey_juvenile)]
        self._predator_juvenile = [value * self.juvenile_predator_fraction for value in self._predator]
        self._predator_adult = [value - juvenile for value, juvenile in zip(self._predator, self._predator_juvenile)]

        self._history: List[Dict[str, Any]] = []
        self._input_overrides: Dict[str, BioSignal] = {}
        self._outputs: Dict[str, BioSignal] = {}
        self._publish_outputs()

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "prey_initial_population": self._scalar_input("count", "Initial prey population count."),
            "predator_initial_population": self._scalar_input("count", "Initial predator population count."),
            "prey_growth_rate": self._scalar_input("1/day", "Intrinsic prey growth rate."),
            "prey_carrying_capacity": self._scalar_input("count", "Prey carrying capacity."),
            "attack_rate": self._scalar_input("1/(count*day)", "Holling type II attack rate."),
            "handling_time": self._scalar_input("day", "Handling time per prey item."),
            "predator_mortality_rate": self._scalar_input("1/day", "Predator mortality rate."),
            "predator_conversion_efficiency": self._scalar_input("predator/prey", "Predator gain per prey consumed."),
            "food_resource_index": self._scalar_input("dimensionless", "Relative food/resource multiplier."),
            "seasonal_forcing": SignalSpec.record(
                schema={"amplitude": "float", "period": "float", "phase": "float"},
                accepted_profiles=(
                    AcceptedSignalProfile(
                        signal_type="record",
                        schema={"amplitude": "float", "period": "float", "phase": "float"},
                    ),
                ),
                description="Seasonal forcing settings for growth, capacity, and mortality.",
            ),
            "disease_parameters": SignalSpec.record(
                schema={"transmission_rate": "float", "recovery_rate": "float", "mortality_rate": "float"},
                accepted_profiles=(
                    AcceptedSignalProfile(
                        signal_type="record",
                        schema={"transmission_rate": "float", "recovery_rate": "float", "mortality_rate": "float"},
                    ),
                ),
                description="Disease rates for susceptible/infected compartments.",
            ),
            "migration_matrix": SignalSpec.record(
                schema={"matrix": "json"},
                accepted_profiles=(AcceptedSignalProfile(signal_type="record", schema={"matrix": "json"}),),
                description="Patch migration matrix or migration-rate configuration.",
            ),
            "patch_configuration": SignalSpec.record(
                schema={"patches": "json"},
                accepted_profiles=(AcceptedSignalProfile(signal_type="record", schema={"patches": "json"}),),
                description="Patch-specific carrying capacities, food resources, and initial populations.",
            ),
        }

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            "visualisation_payload": SignalSpec.record(
                schema={"payload": "json"},
                description="Internal history payload for the sibling visualisation model.",
            ),
            "prey_population_state": SignalSpec.record(
                schema={"role": "str", "label": "str", "count": "float", "patches": "json", "t": "float"},
                emitted_unit="count",
                description="Current prey population state.",
            ),
            "predator_population_state": SignalSpec.record(
                schema={"role": "str", "label": "str", "count": "float", "patches": "json", "t": "float"},
                emitted_unit="count",
                description="Current predator population state.",
            ),
            "population_timeseries": SignalSpec.record(
                schema={"points": "json", "time_unit": "str"},
                emitted_unit="count",
                description="Population trajectory with prey and predator counts.",
            ),
            "equilibrium_summary": SignalSpec.record(
                schema={"prey_equilibrium": "float", "predator_equilibrium": "float", "method": "str"},
                emitted_unit="count",
                description="Rosenzweig-MacArthur coexistence equilibrium summary.",
            ),
            "stability_summary": SignalSpec.record(
                schema={"classification": "str", "prey_growth_margin": "float", "predator_replacement_margin": "float"},
                emitted_unit="dimensionless",
                description="Local qualitative stability indicators.",
            ),
            "extinction_risk": SignalSpec.record(
                schema={"prey": "float", "predator": "float", "joint": "float"},
                emitted_unit="fraction",
                description="Threshold-based extinction risk scores.",
            ),
            "threshold_crossings": SignalSpec.record(
                schema={"prey_below_threshold": "bool", "predator_below_threshold": "bool", "events": "json"},
                emitted_unit="dimensionless",
                description="Recent threshold crossing events.",
            ),
            "scenario_summary": SignalSpec.record(
                schema={"mechanisms": "json", "units": "json", "labels": "json"},
                description="Enabled mechanisms, public units, and display labels.",
            ),
        }

    @staticmethod
    def _scalar_input(unit: str, description: str) -> SignalSpec:
        return SignalSpec.scalar(
            dtype="float64",
            accepted_profiles=(
                AcceptedSignalProfile(
                    signal_type="scalar",
                    dtype="float64",
                    accepted_units=(unit,),
                    description=description,
                ),
            ),
            description=description,
        )

    @staticmethod
    def _fraction(value: float, upper: float = 1.0) -> float:
        return max(0.0, min(upper, float(value)))

    @staticmethod
    def _expand_series(values: Sequence[float] | None, count: int, default: float) -> List[float]:
        if values is None:
            return [max(0.0, float(default)) for _ in range(count)]
        expanded = [max(0.0, float(value)) for value in values]
        if not expanded:
            expanded = [max(0.0, float(default))]
        while len(expanded) < count:
            expanded.append(expanded[-1])
        return expanded[:count]

    @staticmethod
    def _normalize_matrix(matrix: Sequence[Sequence[float]] | None, count: int) -> List[List[float]] | None:
        if matrix is None:
            return None
        normalized: List[List[float]] = []
        for row in matrix[:count]:
            values = [max(0.0, float(value)) for value in row[:count]]
            while len(values) < count:
                values.append(0.0)
            normalized.append(values)
        while len(normalized) < count:
            normalized.append([0.0 for _ in range(count)])
        return normalized

    def set_inputs(self, inputs: dict[str, BioSignal]) -> None:
        self._input_overrides = dict(inputs or {})
        self._apply_input_overrides(reset_initial_state=self._time <= 0.0 and not self._history)

    def _number_input(self, name: str) -> float | None:
        signal = self._input_overrides.get(name)
        if signal is None:
            return None
        return coerce_float(signal)

    def _record_input(self, name: str) -> dict[str, Any]:
        signal = self._input_overrides.get(name)
        if signal is None:
            return {}
        value = unwrap_payload(signal)
        return value if isinstance(value, dict) else {}

    def _apply_input_overrides(self, *, reset_initial_state: bool) -> None:
        for field in (
            "prey_growth_rate",
            "prey_carrying_capacity",
            "attack_rate",
            "handling_time",
            "predator_mortality_rate",
            "predator_conversion_efficiency",
            "food_resource_index",
        ):
            value = self._number_input(field)
            if value is not None:
                setattr(self, field, max(0.0, value))
        self.prey_carrying_capacity = max(EPS, self.prey_carrying_capacity)

        prey_initial = self._number_input("prey_initial_population")
        predator_initial = self._number_input("predator_initial_population")
        if prey_initial is not None:
            self.prey_initial = [max(0.0, prey_initial) for _ in range(self.patch_count)]
            if reset_initial_state:
                self._prey = list(self.prey_initial)
        if predator_initial is not None:
            self.predator_initial = [max(0.0, predator_initial) for _ in range(self.patch_count)]
            if reset_initial_state:
                self._predator = list(self.predator_initial)

        seasonal = self._record_input("seasonal_forcing")
        if seasonal:
            self.seasonal_amplitude = self._fraction(seasonal.get("amplitude", self.seasonal_amplitude), upper=0.95)
            self.seasonal_period = max(EPS, float(seasonal.get("period", self.seasonal_period)))
            self.seasonal_phase = float(seasonal.get("phase", self.seasonal_phase))

        disease = self._record_input("disease_parameters")
        if disease:
            self.disease_enabled = True
            self.prey_transmission_rate = max(0.0, float(disease.get("prey_transmission_rate", disease.get("transmission_rate", self.prey_transmission_rate))))
            self.predator_transmission_rate = max(0.0, float(disease.get("predator_transmission_rate", disease.get("transmission_rate", self.predator_transmission_rate))))
            self.prey_recovery_rate = max(0.0, float(disease.get("prey_recovery_rate", disease.get("recovery_rate", self.prey_recovery_rate))))
            self.predator_recovery_rate = max(0.0, float(disease.get("predator_recovery_rate", disease.get("recovery_rate", self.predator_recovery_rate))))
            self.prey_disease_mortality_rate = max(0.0, float(disease.get("prey_mortality_rate", disease.get("mortality_rate", self.prey_disease_mortality_rate))))
            self.predator_disease_mortality_rate = max(0.0, float(disease.get("predator_mortality_rate", disease.get("mortality_rate", self.predator_disease_mortality_rate))))

        migration = self._record_input("migration_matrix")
        if migration:
            self.migration_rate = max(0.0, float(migration.get("migration_rate", self.migration_rate)))
            if isinstance(migration.get("matrix"), list):
                self.migration_matrix = self._normalize_matrix(migration["matrix"], self.patch_count)

        patch_config = self._record_input("patch_configuration")
        patches = patch_config.get("patches") if patch_config else None
        if isinstance(patches, dict):
            self.patch_food_resources = self._expand_series(patches.get("food_resource_index"), self.patch_count, 1.0)
            self.patch_carrying_capacities = self._expand_series(
                patches.get("prey_carrying_capacity"), self.patch_count, self.prey_carrying_capacity
            )
            if reset_initial_state:
                if "prey_initial_population" in patches:
                    self._prey = self._expand_series(patches["prey_initial_population"], self.patch_count, 0.0)
                if "predator_initial_population" in patches:
                    self._predator = self._expand_series(patches["predator_initial_population"], self.patch_count, 0.0)

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.reset()

    def reset(self) -> None:
        self._time = 0.0
        self._rng = random.Random(self.random_seed)
        self._prey = [max(0.0, float(value)) for value in self.prey_initial]
        self._predator = [max(0.0, float(value)) for value in self.predator_initial]
        self._history = []
        self._publish_outputs()

    def get_state(self) -> dict[str, Any]:
        return {
            "time": self._time,
            "prey": list(self._prey),
            "predator": list(self._predator),
            "prey_total": sum(self._prey),
            "predator_total": sum(self._predator),
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def _seasonal_signal(self, t: float) -> float:
        return math.sin(2.0 * math.pi * (t + self.seasonal_phase) / self.seasonal_period)

    def _effective_carrying_capacity(self, t: float, patch_idx: int) -> float:
        seasonal = self._seasonal_signal(t)
        resource = max(0.0, self.food_resource_index * self.patch_food_resources[patch_idx])
        capacity_modifier = max(0.05, 1.0 + 0.5 * self.seasonal_amplitude * seasonal)
        return max(EPS, self.patch_carrying_capacities[patch_idx] * resource * capacity_modifier)

    def _environment(self, t: float, patch_idx: int) -> tuple[float, float, float]:
        seasonal = self._seasonal_signal(t)
        growth = self.prey_growth_rate * max(0.05, 1.0 + self.seasonal_amplitude * seasonal)
        carrying_capacity = self._effective_carrying_capacity(t, patch_idx)
        mortality = self.predator_mortality_rate * max(0.05, 1.0 - 0.25 * self.seasonal_amplitude * seasonal)
        return growth, carrying_capacity, mortality

    def _predation(self, prey: float, predator: float) -> float:
        if prey <= 0.0 or predator <= 0.0 or self.attack_rate <= 0.0:
            return 0.0
        return self.attack_rate * prey * predator / (1.0 + self.attack_rate * self.handling_time * prey)

    def _disease_losses(self, patch_idx: int) -> tuple[float, float]:
        if not self.disease_enabled:
            return 0.0, 0.0
        return (
            self.prey_disease_mortality_rate * self._prey_infected[patch_idx],
            self.predator_disease_mortality_rate * self._predator_infected[patch_idx],
        )

    def _migration_delta(self, values: Sequence[float]) -> List[float]:
        if self.patch_count <= 1:
            return [0.0 for _ in values]
        if self.migration_matrix is not None:
            deltas = []
            for idx, value in enumerate(values):
                outgoing = sum(self.migration_matrix[idx][j] for j in range(self.patch_count) if j != idx) * value
                incoming = sum(self.migration_matrix[j][idx] * values[j] for j in range(self.patch_count) if j != idx)
                deltas.append(incoming - outgoing)
            return deltas
        if self.migration_rate <= 0.0:
            return [0.0 for _ in values]
        deltas = []
        for idx, value in enumerate(values):
            others = [values[j] for j in range(self.patch_count) if j != idx]
            mean_other = sum(others) / max(1, len(others))
            deltas.append(self.migration_rate * (mean_other - value))
        return deltas

    def _derivatives(self, prey: Sequence[float], predator: Sequence[float], t: float) -> tuple[List[float], List[float], List[float]]:
        prey_migration = self._migration_delta(prey)
        predator_migration = self._migration_delta(predator)
        dprey: List[float] = []
        dpredator: List[float] = []
        predation_rates: List[float] = []
        for idx, (prey_count, predator_count) in enumerate(zip(prey, predator)):
            growth_rate, carrying_capacity, mortality_rate = self._environment(t, idx)
            predation = self._predation(prey_count, predator_count)
            prey_disease_loss, predator_disease_loss = self._disease_losses(idx)
            growth = growth_rate * prey_count * (1.0 - prey_count / carrying_capacity)
            d_n = growth - predation - prey_disease_loss + prey_migration[idx]
            d_p = (
                self.predator_conversion_efficiency * predation
                - mortality_rate * predator_count
                - predator_disease_loss
                + predator_migration[idx]
            )
            if self.noise_strength > 0.0:
                d_n += self._rng.gauss(0.0, self.noise_strength) * max(1.0, prey_count)
                d_p += self._rng.gauss(0.0, self.noise_strength) * max(1.0, predator_count)
            dprey.append(d_n)
            dpredator.append(d_p)
            predation_rates.append(predation)
        return dprey, dpredator, predation_rates

    def advance_window(
        self,
        start: float | None = None,
        end: float | None = None,
        inputs: dict[str, BioSignal] | None = None,
    ) -> dict[str, BioSignal]:
        if inputs:
            self.set_inputs(inputs)
        else:
            self._apply_input_overrides(reset_initial_state=False)

        if end is None:
            end = self._time + float(getattr(self, "communication_step", self.integration_step) or self.integration_step)
        target = float(end)
        if target <= self._time:
            return dict(self._outputs)

        while self._time < target - 1e-12:
            dt = min(self.integration_step, target - self._time)
            dprey, dpredator, predation_rates = self._derivatives(self._prey, self._predator, self._time)
            self._prey = [max(0.0, value + dt * delta) for value, delta in zip(self._prey, dprey)]
            self._predator = [max(0.0, value + dt * delta) for value, delta in zip(self._predator, dpredator)]
            self._update_disease(dt)
            self._update_stage_structure(dt, predation_rates)
            self._time += dt

        self._append_history()
        self._publish_outputs()
        return dict(self._outputs)

    def _update_disease(self, dt: float) -> None:
        if not self.disease_enabled:
            self._prey_susceptible = list(self._prey)
            self._predator_susceptible = list(self._predator)
            self._prey_infected = [0.0 for _ in self._prey]
            self._predator_infected = [0.0 for _ in self._predator]
            return
        for idx in range(self.patch_count):
            prey_total = max(EPS, self._prey[idx])
            predator_total = max(EPS, self._predator[idx])
            prey_infections = self.prey_transmission_rate * self._prey_susceptible[idx] * self._prey_infected[idx] / prey_total
            predator_infections = (
                self.predator_transmission_rate
                * self._predator_susceptible[idx]
                * self._predator_infected[idx]
                / predator_total
            )
            prey_recoveries = self.prey_recovery_rate * self._prey_infected[idx]
            predator_recoveries = self.predator_recovery_rate * self._predator_infected[idx]
            prey_deaths = self.prey_disease_mortality_rate * self._prey_infected[idx]
            predator_deaths = self.predator_disease_mortality_rate * self._predator_infected[idx]
            self._prey_infected[idx] = max(
                0.0, self._prey_infected[idx] + dt * (prey_infections - prey_recoveries - prey_deaths)
            )
            self._predator_infected[idx] = max(
                0.0,
                self._predator_infected[idx] + dt * (predator_infections - predator_recoveries - predator_deaths),
            )
            self._prey_infected[idx] = min(self._prey_infected[idx], self._prey[idx])
            self._predator_infected[idx] = min(self._predator_infected[idx], self._predator[idx])
            self._prey_susceptible[idx] = max(0.0, self._prey[idx] - self._prey_infected[idx])
            self._predator_susceptible[idx] = max(0.0, self._predator[idx] - self._predator_infected[idx])

    def _update_stage_structure(self, dt: float, predation_rates: Sequence[float]) -> None:
        if not self.age_structure_enabled:
            self._prey_juvenile = [value * self.juvenile_prey_fraction for value in self._prey]
            self._prey_adult = [value - juvenile for value, juvenile in zip(self._prey, self._prey_juvenile)]
            self._predator_juvenile = [value * self.juvenile_predator_fraction for value in self._predator]
            self._predator_adult = [value - juvenile for value, juvenile in zip(self._predator, self._predator_juvenile)]
            return
        for idx in range(self.patch_count):
            prey_births = max(0.0, self.prey_growth_rate * self._prey[idx] * dt)
            predator_births = max(0.0, self.predator_conversion_efficiency * predation_rates[idx] * dt)
            prey_maturation = min(self._prey_juvenile[idx], self.prey_maturation_rate * self._prey_juvenile[idx] * dt)
            predator_maturation = min(
                self._predator_juvenile[idx], self.predator_maturation_rate * self._predator_juvenile[idx] * dt
            )
            self._prey_juvenile[idx] = max(
                0.0,
                self._prey_juvenile[idx]
                + prey_births
                - prey_maturation
                - self.juvenile_prey_mortality_rate * self._prey_juvenile[idx] * dt,
            )
            self._prey_adult[idx] = max(0.0, self._prey_adult[idx] + prey_maturation)
            self._predator_juvenile[idx] = max(
                0.0,
                self._predator_juvenile[idx]
                + predator_births
                - predator_maturation
                - self.juvenile_predator_mortality_rate * self._predator_juvenile[idx] * dt,
            )
            self._predator_adult[idx] = max(0.0, self._predator_adult[idx] + predator_maturation)
            self._renormalize_stage(idx, self._prey[idx], self._prey_juvenile, self._prey_adult)
            self._renormalize_stage(idx, self._predator[idx], self._predator_juvenile, self._predator_adult)

    @staticmethod
    def _renormalize_stage(idx: int, total: float, juvenile: List[float], adult: List[float]) -> None:
        stage_total = juvenile[idx] + adult[idx]
        if stage_total <= EPS:
            juvenile[idx] = 0.0
            adult[idx] = max(0.0, total)
            return
        scale = max(0.0, total) / stage_total
        juvenile[idx] *= scale
        adult[idx] *= scale

    def _append_history(self) -> None:
        point = {
            "t": float(self._time),
            "prey": float(sum(self._prey)),
            "predator": float(sum(self._predator)),
            "prey_patches": [float(value) for value in self._prey],
            "predator_patches": [float(value) for value in self._predator],
        }
        self._history.append(point)
        if len(self._history) > HISTORY_POINT_LIMIT:
            self._history = self._history[-HISTORY_POINT_LIMIT:]

    def _equilibrium_summary(self) -> dict[str, Any]:
        replacement = self.predator_conversion_efficiency - self.predator_mortality_rate * self.handling_time
        if self.attack_rate <= EPS or replacement <= EPS:
            return {"prey_equilibrium": -1.0, "predator_equilibrium": 0.0, "method": "no_positive_coexistence"}
        prey_equilibrium = self.predator_mortality_rate / (self.attack_rate * replacement)
        carrying_capacity = self.prey_carrying_capacity * max(EPS, self.food_resource_index)
        predation_per_predator = self.attack_rate * prey_equilibrium / (
            1.0 + self.attack_rate * self.handling_time * prey_equilibrium
        )
        predator_equilibrium = max(
            0.0,
            self.prey_growth_rate
            * (1.0 - prey_equilibrium / max(EPS, carrying_capacity))
            * prey_equilibrium
            / max(EPS, predation_per_predator),
        )
        return {
            "prey_equilibrium": float(prey_equilibrium),
            "predator_equilibrium": float(predator_equilibrium),
            "method": "holling_type_ii_coexistence",
        }

    def _stability_summary(self) -> dict[str, Any]:
        prey_total = sum(self._prey)
        carrying_capacity = self.prey_carrying_capacity * max(EPS, self.food_resource_index)
        prey_growth_margin = carrying_capacity - prey_total
        replacement = (
            self.predator_conversion_efficiency
            * self.attack_rate
            * max(0.0, prey_total)
            / (1.0 + self.attack_rate * self.handling_time * max(0.0, prey_total))
            - self.predator_mortality_rate
        )
        if replacement < 0:
            classification = "predator_decline_likely"
        elif prey_total > carrying_capacity:
            classification = "resource_overshoot"
        else:
            classification = "bounded_coexistence_possible"
        return {
            "classification": classification,
            "prey_growth_margin": float(prey_growth_margin),
            "predator_replacement_margin": float(replacement),
        }

    def _extinction_risk(self) -> dict[str, float]:
        if not self._history:
            prey = sum(self._prey)
            predator = sum(self._predator)
            prey_risk = 1.0 if prey <= self.prey_extinction_threshold else 0.0
            predator_risk = 1.0 if predator <= self.predator_extinction_threshold else 0.0
            return {"prey": prey_risk, "predator": predator_risk, "joint": max(prey_risk, predator_risk)}
        recent = self._history[-min(50, len(self._history)) :]
        prey_hits = sum(1 for point in recent if point["prey"] <= self.prey_extinction_threshold)
        predator_hits = sum(1 for point in recent if point["predator"] <= self.predator_extinction_threshold)
        denom = float(len(recent))
        prey_risk = prey_hits / denom
        predator_risk = predator_hits / denom
        return {"prey": prey_risk, "predator": predator_risk, "joint": max(prey_risk, predator_risk)}

    def _threshold_crossings(self) -> dict[str, Any]:
        events = []
        if self._history:
            previous = self._history[-2] if len(self._history) >= 2 else self._history[-1]
            current = self._history[-1]
            if previous["prey"] > self.prey_extinction_threshold >= current["prey"]:
                events.append({"t": current["t"], "target": "prey", "threshold": self.prey_extinction_threshold})
            if previous["predator"] > self.predator_extinction_threshold >= current["predator"]:
                events.append({"t": current["t"], "target": "predator", "threshold": self.predator_extinction_threshold})
        return {
            "prey_below_threshold": bool(sum(self._prey) <= self.prey_extinction_threshold),
            "predator_below_threshold": bool(sum(self._predator) <= self.predator_extinction_threshold),
            "events": events,
        }

    def _scenario_summary(self) -> dict[str, Any]:
        return {
            "mechanisms": {
                "carrying_capacity": True,
                "holling_type_ii_predation": True,
                "seasonality": self.seasonal_amplitude > 0.0,
                "noise": self.noise_strength > 0.0,
                "food_limitation": self.food_resource_index < 1.0 or any(value != 1.0 for value in self.patch_food_resources),
                "age_structure": self.age_structure_enabled,
                "disease": self.disease_enabled,
                "migration": self.patch_count > 1 and (self.migration_rate > 0.0 or self.migration_matrix is not None),
                "spatial_patches": self.patch_count,
            },
            "units": {
                "population": "count",
                "time": self.time_unit,
                "prey_growth_rate": "1/day",
                "prey_carrying_capacity": "count",
                "attack_rate": "1/(count*day)",
                "handling_time": "day",
                "predator_mortality_rate": "1/day",
                "predator_conversion_efficiency": "predator/prey",
                "food_resource_index": "dimensionless",
            },
            "labels": {"prey": self.prey_label, "predator": self.predator_label},
        }

    def _publish_outputs(self) -> None:
        source = getattr(self, "_world_name", self.__class__.__name__)
        specs = self.outputs()
        prey_total = float(sum(self._prey))
        predator_total = float(sum(self._predator))
        self._outputs = {
            "visualisation_payload": _make_signal(
                source=source,
                name="visualisation_payload",
                value={"payload": self._visualisation_payload()},
                emitted_at=float(self._time),
                spec=specs["visualisation_payload"],
            ),
            "prey_population_state": _make_signal(
                source=source,
                name="prey_population_state",
                value={
                    "role": "prey",
                    "label": self.prey_label,
                    "count": prey_total,
                    "patches": [float(value) for value in self._prey],
                    "t": float(self._time),
                },
                emitted_at=float(self._time),
                spec=specs["prey_population_state"],
            ),
            "predator_population_state": _make_signal(
                source=source,
                name="predator_population_state",
                value={
                    "role": "predator",
                    "label": self.predator_label,
                    "count": predator_total,
                    "patches": [float(value) for value in self._predator],
                    "t": float(self._time),
                },
                emitted_at=float(self._time),
                spec=specs["predator_population_state"],
            ),
            "population_timeseries": _make_signal(
                source=source,
                name="population_timeseries",
                value={"points": list(self._history), "time_unit": self.time_unit},
                emitted_at=float(self._time),
                spec=specs["population_timeseries"],
            ),
            "equilibrium_summary": _make_signal(
                source=source,
                name="equilibrium_summary",
                value=self._equilibrium_summary(),
                emitted_at=float(self._time),
                spec=specs["equilibrium_summary"],
            ),
            "stability_summary": _make_signal(
                source=source,
                name="stability_summary",
                value=self._stability_summary(),
                emitted_at=float(self._time),
                spec=specs["stability_summary"],
            ),
            "extinction_risk": _make_signal(
                source=source,
                name="extinction_risk",
                value=self._extinction_risk(),
                emitted_at=float(self._time),
                spec=specs["extinction_risk"],
            ),
            "threshold_crossings": _make_signal(
                source=source,
                name="threshold_crossings",
                value=self._threshold_crossings(),
                emitted_at=float(self._time),
                spec=specs["threshold_crossings"],
            ),
            "scenario_summary": _make_signal(
                source=source,
                name="scenario_summary",
                value=self._scenario_summary(),
                emitted_at=float(self._time),
                spec=specs["scenario_summary"],
            ),
        }

    def _visualisation_payload(self) -> Dict[str, Any]:
        return {
            "time_unit": self.time_unit,
            "prey_label": self.prey_label,
            "predator_label": self.predator_label,
            "equilibrium_summary": self._equilibrium_summary(),
            "stability_summary": self._stability_summary(),
            "extinction_risk": self._extinction_risk(),
            "history": list(self._history),
        }

    def visualize(self) -> Optional["VisualSpec" | List["VisualSpec"]]:
        return None

    def _timeseries_visual(self) -> "VisualSpec":
        points = self._visual_history()
        return {
            "render": "timeseries",
            "description": "Prey and predator population counts under the configured Rosenzweig-MacArthur scenario.",
            "data": {
                "title": "Population Trajectories",
                "x_unit": self.time_unit,
                "y_unit": "count",
                "series": [
                    {"name": self.prey_label, "points": [[point["t"], point["prey"]] for point in points]},
                    {"name": self.predator_label, "points": [[point["t"], point["predator"]] for point in points]},
                ],
            },
        }

    def _visual_history(self, limit: int = VISUAL_POINT_LIMIT) -> List[Dict[str, Any]]:
        if len(self._history) <= limit:
            return list(self._history)
        if limit <= 1:
            return [self._history[-1]]

        max_index = len(self._history) - 1
        sampled: List[Dict[str, Any]] = []
        last_index = -1
        for index in range(limit):
            history_index = round(index * max_index / (limit - 1))
            if history_index != last_index:
                sampled.append(self._history[history_index])
                last_index = history_index
        return sampled

    def _phase_visual(self) -> "VisualSpec":
        svg = self._phase_svg()
        encoded = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        return {
            "render": "image",
            "description": "Phase portrait with trajectory, effective carrying capacity, and coexistence equilibrium marker.",
            "data": {
                "title": "Rosenzweig-MacArthur Phase Portrait",
                "src": f"data:image/svg+xml;base64,{encoded}",
                "alt": "Rosenzweig-MacArthur phase portrait",
                "width": 720,
                "height": 420,
            },
        }

    def _functional_response_visual(self) -> "VisualSpec":
        svg = self._functional_response_svg()
        encoded = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        return {
            "render": "image",
            "description": "Holling type II functional response showing predation saturation and the current prey density.",
            "data": {
                "title": "Holling Type II Functional Response",
                "src": f"data:image/svg+xml;base64,{encoded}",
                "alt": "Holling type II functional response curve",
                "width": 720,
                "height": 360,
            },
        }

    def _ecology_state_table_visual(self) -> "VisualSpec":
        risk = self._extinction_risk()
        rows = []
        for idx in range(self.patch_count):
            rows.append(
                [
                    idx + 1,
                    f"{self._prey[idx]:.6g}",
                    f"{self._predator[idx]:.6g}",
                    f"{self._effective_carrying_capacity(self._time, idx):.6g}",
                    f"{self.patch_food_resources[idx]:.6g}",
                    f"{self._prey_infected[idx]:.6g}",
                    f"{self._predator_infected[idx]:.6g}",
                ]
            )
        rows.extend(
            [
                ["risk", f"prey={risk['prey']:.3g}", f"predator={risk['predator']:.3g}", f"joint={risk['joint']:.3g}", "", "", ""],
                ["mechanisms", self._mechanism_label(), "", "", "", "", ""],
            ]
        )
        return {
            "render": "table",
            "description": "Patch resource state, disease compartments, and threshold-based extinction risk.",
            "data": {
                "title": "Ecology State and Risk",
                "columns": [
                    "Patch",
                    "Prey count",
                    "Predator count",
                    "Effective K",
                    "Food index",
                    "Prey infected",
                    "Predator infected",
                ],
                "rows": rows,
            },
        }

    def _mechanism_label(self) -> str:
        mechanisms = self._scenario_summary()["mechanisms"]
        return ", ".join(name for name, enabled in mechanisms.items() if enabled is True)

    def _phase_svg(self) -> str:
        width, height = 720, 420
        left, right, top, bottom = 64, 24, 32, 50
        plot_w = width - left - right
        plot_h = height - top - bottom
        visual_history = self._visual_history()
        points = [(p["prey"], p["predator"]) for p in visual_history]
        if len(points) < 2:
            points = [(sum(self._prey), sum(self._predator)), (sum(self._prey), sum(self._predator))]
        equilibrium = self._equilibrium_summary()
        carrying_capacity = max(self._effective_carrying_capacity(self._time, i) for i in range(self.patch_count))
        eq_prey = max(0.0, equilibrium["prey_equilibrium"])
        eq_predator = max(0.0, equilibrium["predator_equilibrium"])
        max_x = max([x for x, _ in points] + [carrying_capacity, eq_prey, 1.0]) * 1.08
        max_y = max([y for _, y in points] + [eq_predator, 1.0]) * 1.08

        def sx(value: float) -> float:
            return left + plot_w * value / max(EPS, max_x)

        def sy(value: float) -> float:
            return top + plot_h * (1.0 - value / max(EPS, max_y))

        path = " ".join(("M" if i == 0 else "L") + f" {sx(x):.2f},{sy(y):.2f}" for i, (x, y) in enumerate(points))
        k_x = sx(carrying_capacity)
        start_x, start_y = points[0]
        end_x, end_y = points[-1]
        eq = "" if equilibrium["method"] == "no_positive_coexistence" else f'<circle cx="{sx(eq_prey):.2f}" cy="{sy(eq_predator):.2f}" r="5" fill="#facc15"/><text x="{sx(eq_prey)+8:.2f}" y="{sy(eq_predator)-8:.2f}" fill="#e5e7eb" font-size="12">coexistence equilibrium</text>'
        return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">
  <rect width="{width}" height="{height}" fill="#111827" rx="10"/>
  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#64748b"/>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#64748b"/>
  <line x1="{k_x:.2f}" y1="{top}" x2="{k_x:.2f}" y2="{height-bottom}" stroke="#22c55e" stroke-dasharray="6 5"/>
  <path d="{path}" fill="none" stroke="#38bdf8" stroke-width="3"/>
  <circle cx="{sx(start_x):.2f}" cy="{sy(start_y):.2f}" r="4" fill="#22c55e"/>
  <circle cx="{sx(end_x):.2f}" cy="{sy(end_y):.2f}" r="4" fill="#ef4444"/>
  {eq}
  <text x="{width/2}" y="23" text-anchor="middle" fill="#f8fafc" font-size="16" font-family="sans-serif">Rosenzweig-MacArthur Phase Portrait</text>
  <text x="{k_x + 6:.2f}" y="{top + 16}" fill="#86efac" font-size="12">effective K</text>
  <text x="{left + 8}" y="{top + 16}" fill="#86efac" font-size="12">start</text>
  <text x="{left + 8}" y="{top + 32}" fill="#fca5a5" font-size="12">end</text>
  <text x="{width/2}" y="{height-14}" text-anchor="middle" fill="#cbd5e1" font-size="13">Prey population (count)</text>
  <text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" fill="#cbd5e1" font-size="13">Predator population (count)</text>
</svg>"""

    def _functional_response_svg(self) -> str:
        width, height = 720, 360
        left, right, top, bottom = 64, 24, 30, 48
        plot_w = width - left - right
        plot_h = height - top - bottom
        max_prey = max(max(sum(self._prey), self.prey_carrying_capacity) * 1.25, 1.0)
        curve = []
        max_rate = 0.0
        for i in range(80):
            prey = max_prey * i / 79
            rate = self._predation(prey, 1.0)
            curve.append((prey, rate))
            max_rate = max(max_rate, rate)
        max_rate = max(max_rate, EPS)

        def sx(value: float) -> float:
            return left + plot_w * value / max_prey

        def sy(value: float) -> float:
            return top + plot_h * (1.0 - value / max_rate)

        path = " ".join(("M" if i == 0 else "L") + f" {sx(x):.2f},{sy(y):.2f}" for i, (x, y) in enumerate(curve))
        current_x = sx(sum(self._prey))
        return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">
  <rect width="{width}" height="{height}" fill="#111827" rx="10"/>
  <line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#64748b"/>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#64748b"/>
  <path d="{path}" fill="none" stroke="#f97316" stroke-width="3"/>
  <line x1="{current_x:.2f}" y1="{top}" x2="{current_x:.2f}" y2="{height-bottom}" stroke="#38bdf8" stroke-dasharray="5 5"/>
  <text x="{width/2}" y="22" text-anchor="middle" fill="#f8fafc" font-size="16" font-family="sans-serif">Holling Type II Predation Saturation</text>
  <text x="{current_x + 6:.2f}" y="{top + 18}" fill="#bae6fd" font-size="12">current prey</text>
  <text x="{width/2}" y="{height-14}" text-anchor="middle" fill="#cbd5e1" font-size="13">Prey density (count)</text>
  <text x="18" y="{height/2}" transform="rotate(-90 18 {height/2})" text-anchor="middle" fill="#cbd5e1" font-size="13">Prey consumed per predator per day</text>
</svg>"""

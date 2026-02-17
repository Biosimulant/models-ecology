# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Predator-prey interaction using Lotka-Volterra-style dynamics."""
from __future__ import annotations

import random
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from bsim import BioWorld
    from bsim.visuals import VisualSpec

from bsim import BioModule
from bsim.signals import BioSignal, SignalMetadata


class PredatorPreyInteraction(BioModule):
    """Models predator-prey interactions using Lotka-Volterra-style dynamics.

    Receives population states from predator and prey modules, computes kills,
    and emits predation signals to prey and food signals to predators.

    Parameters:
        predation_rate: Base rate of successful hunts (predator * prey * rate = encounters).
        conversion_efficiency: Fraction of prey biomass converted to predator food.
        satiation_factor: Predators hunt less when well-fed (0 = no effect, 1 = strong effect).
        min_prey_for_hunt: Minimum prey count before hunting is possible.
        seed: Random seed for reproducibility.
    """

    def __init__(
        self,
        predation_rate: float = 0.001,
        conversion_efficiency: float = 0.1,
        satiation_factor: float = 0.0,
        min_prey_for_hunt: int = 0,
        seed: Optional[int] = None,
        min_dt: float = 1.0,
    ) -> None:
        self.min_dt = min_dt
        self.predation_rate = predation_rate
        self.conversion_efficiency = conversion_efficiency
        self.satiation_factor = satiation_factor
        self.min_prey_for_hunt = min_prey_for_hunt
        self.seed = seed
        self._rng = random.Random(seed)

        self._prey_count: int = 0
        self._prey_species: str = "Prey"
        self._predator_count: int = 0
        self._predator_species: str = "Predator"
        self._time: float = 0.0
        self._history: List[Dict[str, Any]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"prey_state", "predator_state"}

    def outputs(self) -> Set[str]:
        return {"predation", "food_gained"}

    def reset(self) -> None:
        """Reset interaction state."""
        self._rng = random.Random(self.seed)
        self._prey_count = 0
        self._predator_count = 0
        self._time = 0.0
        self._history = []
        self._outputs = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        prey = signals.get("prey_state")
        if prey is not None and isinstance(prey.value, dict):
            self._prey_count = int(prey.value.get("count", 0))
            self._prey_species = str(prey.value.get("species", "Prey"))
        predator = signals.get("predator_state")
        if predator is not None and isinstance(predator.value, dict):
            self._predator_count = int(predator.value.get("count", 0))
            self._predator_species = str(predator.value.get("species", "Predator"))

    def advance_to(self, t: float) -> None:
        dt = t - self._time if t > self._time else self.min_dt
        self._time = t

        kills = 0
        food_gained = 0.0

        if (self._prey_count >= self.min_prey_for_hunt and
            self._predator_count > 0 and
            self._prey_count > 0):

            # Lotka-Volterra style: encounters = predators * prey * rate
            expected_kills = (
                self.predation_rate *
                self._predator_count *
                self._prey_count *
                dt
            )

            # Apply satiation (fewer kills when predators are well-fed)
            if self.satiation_factor > 0:
                ratio = self._predator_count / max(1, self._prey_count)
                satiation_mult = max(0.1, 1 - self.satiation_factor * ratio * 10)
                expected_kills *= satiation_mult

            # Stochastic kills
            if expected_kills > 0:
                # Poisson-like sampling
                if expected_kills > 20:
                    kills = max(0, int(self._rng.gauss(expected_kills, expected_kills ** 0.5)))
                else:
                    import math
                    L = math.exp(-expected_kills)
                    k = 0
                    p = 1.0
                    while p > L:
                        k += 1
                        p *= self._rng.random()
                    kills = k - 1

            # Can't kill more than exist
            kills = min(kills, self._prey_count)

            # Food gained by predators
            food_gained = kills * self.conversion_efficiency

        # Record history
        self._history.append({
            "t": t,
            "kills": kills,
            "food_gained": food_gained,
            "prey_count": self._prey_count,
            "predator_count": self._predator_count,
        })

        source_name = getattr(self, "_world_name", self.__class__.__name__)
        outputs: Dict[str, BioSignal] = {}
        if kills > 0:
            outputs["predation"] = BioSignal(
                source=source_name,
                name="predation",
                value={
                    "kills": kills,
                    "predator": self._predator_species,
                    "t": t,
                },
                time=t,
                metadata=SignalMetadata(units=None, description="Predation events", kind="event"),
            )
        if food_gained > 0:
            outputs["food_gained"] = BioSignal(
                source=source_name,
                name="food_gained",
                value=food_gained,
                time=t,
                metadata=SignalMetadata(units=None, description="Food gained", kind="event"),
            )
        self._outputs = outputs

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional["VisualSpec"]:
        """Generate visualization of predation events over time."""
        if not self._history:
            return None

        return {
            "render": "timeseries",
            "data": {
                "series": [
                    {
                        "name": "Kills per Step",
                        "points": [[h["t"], h["kills"]] for h in self._history],
                    },
                    {
                        "name": "Food Gained",
                        "points": [[h["t"], h["food_gained"]] for h in self._history],
                    },
                ],
                "title": f"Predation: {self._predator_species} \u2192 {self._prey_species}",
            },
        }


class CompetitionInteraction(BioModule):
    """Models competition between species for shared resources.

    When multiple species compete for the same food/space, this module
    reduces effective food availability based on competitor populations.

    Parameters:
        competition_coefficient: How strongly competitors affect each other (0-1).
        resource_type: Type of resource competed for ("food", "space", "water").
    """

    def __init__(
        self,
        competition_coefficient: float = 0.5,
        resource_type: str = "food",
        min_dt: float = 1.0,
    ) -> None:
        self.min_dt = min_dt
        self.competition_coefficient = competition_coefficient
        self.resource_type = resource_type

        self._populations: Dict[str, int] = {}  # species -> count
        self._time: float = 0.0
        self._history: List[Dict[str, Any]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"population_state"}

    def outputs(self) -> Set[str]:
        return {"competition"}

    def reset(self) -> None:
        """Reset competition state."""
        self._populations = {}
        self._time = 0.0
        self._history = []
        self._outputs = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        signal = signals.get("population_state")
        if signal is None or not isinstance(signal.value, dict):
            return
        species = str(signal.value.get("species", "Unknown"))
        count = int(signal.value.get("count", 0))
        self._populations[species] = count

    def advance_to(self, t: float) -> None:
        self._time = t

        if len(self._populations) < 2:
            return  # Need at least 2 species to compete

        total_pop = sum(self._populations.values())

        # Calculate competition pressure for each species
        entries: List[Dict[str, Any]] = []
        for species, count in self._populations.items():
            if count <= 0:
                continue

            # Competitors = everyone else
            competitors = total_pop - count

            # Competition effect: reduces effective resources
            competition_pressure = (
                competitors * self.competition_coefficient / max(1, total_pop)
            )

            entries.append(
                {
                    "species": species,
                    "pressure": competition_pressure,
                    "resource": self.resource_type,
                    "t": t,
                }
            )

        # Record history
        self._history.append({
            "t": t,
            "total_population": total_pop,
            "populations": dict(self._populations),
        })

        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "competition": BioSignal(
                source=source_name,
                name="competition",
                value=entries,
                time=t,
                metadata=SignalMetadata(units=None, description="Competition pressures", kind="state"),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional["VisualSpec"]:
        """Visualize total competing population over time."""
        if not self._history:
            return None

        return {
            "render": "timeseries",
            "data": {
                "series": [
                    {
                        "name": "Total Competing Pop",
                        "points": [[h["t"], h["total_population"]] for h in self._history],
                    },
                ],
                "title": f"Competition for {self.resource_type}",
            },
        }


class MutualismInteraction(BioModule):
    """Models mutualistic (beneficial) interactions between species.

    Both species benefit from each other's presence, increasing survival
    or reproduction rates.

    Parameters:
        benefit_rate: How much each species benefits from the other (0-1).
        benefit_type: Type of benefit ("food", "protection", "reproduction").
    """

    def __init__(
        self,
        benefit_rate: float = 0.1,
        benefit_type: str = "food",
        min_dt: float = 1.0,
    ) -> None:
        self.min_dt = min_dt
        self.benefit_rate = benefit_rate
        self.benefit_type = benefit_type

        self._species_a_count: int = 0
        self._species_a_name: str = "Species A"
        self._species_b_count: int = 0
        self._species_b_name: str = "Species B"
        self._time: float = 0.0
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"species_a_state", "species_b_state"}

    def outputs(self) -> Set[str]:
        return {"mutualism_benefit"}

    def reset(self) -> None:
        self._species_a_count = 0
        self._species_b_count = 0
        self._time = 0.0
        self._outputs = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        a_state = signals.get("species_a_state")
        if a_state is not None and isinstance(a_state.value, dict):
            self._species_a_count = int(a_state.value.get("count", 0))
            self._species_a_name = str(a_state.value.get("species", "Species A"))
        b_state = signals.get("species_b_state")
        if b_state is not None and isinstance(b_state.value, dict):
            self._species_b_count = int(b_state.value.get("count", 0))
            self._species_b_name = str(b_state.value.get("species", "Species B"))

    def advance_to(self, t: float) -> None:
        self._time = t

        if self._species_a_count > 0 and self._species_b_count > 0:
            # Benefit scales with partner population (with diminishing returns)
            import math

            benefit_to_a = self.benefit_rate * math.log1p(self._species_b_count) / 10
            benefit_to_b = self.benefit_rate * math.log1p(self._species_a_count) / 10
            entries = [
                {
                    "species": self._species_a_name,
                    "benefit": benefit_to_a,
                    "type": self.benefit_type,
                    "t": t,
                },
                {
                    "species": self._species_b_name,
                    "benefit": benefit_to_b,
                    "type": self.benefit_type,
                    "t": t,
                },
            ]
            source_name = getattr(self, "_world_name", self.__class__.__name__)
            self._outputs = {
                "mutualism_benefit": BioSignal(
                    source=source_name,
                    name="mutualism_benefit",
                    value=entries,
                    time=t,
                    metadata=SignalMetadata(units=None, description="Mutualism benefits", kind="event"),
                )
            }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional["VisualSpec"]:
        return None  # Mutualism effects are reflected in population dynamics

# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Stochastic multispecies competition model grounded in Leibovich et al. (2022).

The upstream BioModels record for MODEL2212080001 distributes Python Gillespie
code, not SBML. This Biosimulant package therefore exposes a curated
competition-immigration model with explicit ecological observables instead of a
generic SBML wrapper.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

import biosim
import numpy as np
from biosim.signals import BioSignal, SignalMetadata

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


class Leibovich2022CommunityModel(biosim.BioModule):
    """Competition-immigration model with demographic noise."""

    def __init__(
        self,
        species_count: int = 6,
        carrying_capacity: float = 600.0,
        birth_rate: float = 1.0,
        death_rate: float = 0.35,
        competition_strength: float = 0.7,
        competition_overlap: float = 0.65,
        immigration_rate: float = 0.6,
        initial_abundance: float = 45.0,
        rng_seed: int = 7,
        min_dt: float = 0.05,
    ) -> None:
        if min_dt <= 0:
            raise ValueError("min_dt must be positive")
        if species_count < 2:
            raise ValueError("species_count must be at least 2")

        self.min_dt = float(min_dt)
        self.species_count = int(species_count)
        self.carrying_capacity = float(carrying_capacity)
        self.birth_rate = float(birth_rate)
        self.death_rate = float(death_rate)
        self.competition_strength = float(competition_strength)
        self.competition_overlap = float(competition_overlap)
        self.immigration_rate = float(immigration_rate)
        self.initial_abundance = float(initial_abundance)
        self.rng_seed = int(rng_seed)

        self._time = 0.0
        self._rng = np.random.default_rng(self.rng_seed)
        self._abundances = self._initial_state()
        self._history: List[Dict[str, float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.reset()

    def reset(self) -> None:
        self._time = 0.0
        self._rng = np.random.default_rng(self.rng_seed)
        self._abundances = self._initial_state()
        self._history = []
        self._outputs = {}

    def inputs(self) -> Set[str]:
        return set()

    def outputs(self) -> Set[str]:
        return {"community_state", "diversity_metrics"}

    def advance_to(self, t: float) -> None:
        if t <= self._time:
            return

        current = self._time
        while current < t - 1e-12:
            dt = min(self.min_dt, t - current)
            self._step(dt)
            current += dt
            self._record_state(current)

        self._time = current
        self._publish_outputs(current)

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[List["VisualSpec"]]:
        if not self._history:
            return None
        return [
            self._abundance_visual(),
            self._metrics_visual(),
            self._summary_visual(),
        ]

    def _initial_state(self) -> np.ndarray:
        offsets = np.linspace(1.3, 0.7, self.species_count)
        return np.maximum(0, np.rint(self.initial_abundance * offsets)).astype(int)

    def _step(self, dt: float) -> None:
        total_abundance = float(np.sum(self._abundances))
        next_state = self._abundances.copy()

        immigration_weights = np.linspace(1.0, 1.0 + 0.35 * (self.species_count - 1), self.species_count)
        immigration_weights = immigration_weights / float(np.sum(immigration_weights))

        for idx, abundance in enumerate(self._abundances):
            own_density = float(abundance)
            crowding = (
                (1.0 - self.competition_overlap) * own_density
                + self.competition_overlap * total_abundance
            ) / max(self.carrying_capacity, 1.0)
            births_mean = max(0.0, (self.birth_rate * own_density + self.immigration_rate * immigration_weights[idx]) * dt)
            deaths_mean = max(0.0, (self.death_rate + self.competition_strength * crowding) * own_density * dt)

            births = int(self._rng.poisson(births_mean))
            deaths = int(self._rng.poisson(deaths_mean))
            next_state[idx] = max(0, int(abundance) + births - deaths)

        self._abundances = next_state

    def _record_state(self, t: float) -> None:
        total_abundance = int(np.sum(self._abundances))
        richness = int(np.count_nonzero(self._abundances > 0))
        proportions = self._abundances[self._abundances > 0] / max(float(total_abundance), 1.0)
        shannon = -float(np.sum(proportions * np.log(proportions))) if proportions.size else 0.0
        evenness = shannon / math.log(richness) if richness > 1 else 0.0
        dominant_index = int(np.argmax(self._abundances))

        record: Dict[str, float] = {
            "t": float(t),
            "total_abundance": float(total_abundance),
            "richness": float(richness),
            "shannon_diversity": float(shannon),
            "evenness": float(evenness),
            "dominant_species_index": float(dominant_index + 1),
        }
        for idx, abundance in enumerate(self._abundances, start=1):
            record[f"species_{idx}"] = float(abundance)
        self._history.append(record)

    def _publish_outputs(self, t: float) -> None:
        latest = self._history[-1]
        species_abundances = {
            f"species_{idx}": int(self._abundances[idx - 1])
            for idx in range(1, self.species_count + 1)
        }
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "community_state": BioSignal(
                source=source_name,
                name="community_state",
                value={
                    "species_abundances": species_abundances,
                    "total_abundance": int(latest["total_abundance"]),
                },
                time=float(t),
                metadata=SignalMetadata(
                    units="individuals",
                    description="Species-resolved abundances for the stochastic competition community.",
                    kind="state",
                ),
            ),
            "diversity_metrics": BioSignal(
                source=source_name,
                name="diversity_metrics",
                value={
                    "richness": int(latest["richness"]),
                    "shannon_diversity": float(latest["shannon_diversity"]),
                    "evenness": float(latest["evenness"]),
                    "dominant_species": f"species_{int(latest['dominant_species_index'])}",
                },
                time=float(t),
                metadata=SignalMetadata(
                    units=None,
                    description="Diversity and dominance metrics for the multispecies community.",
                    kind="summary",
                ),
            ),
        }

    def _abundance_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Species-resolved abundance trajectories for the stochastic competition-immigration community.",
            "data": {
                "title": "Community Abundance Trajectories",
                "series": [
                    {
                        "name": f"species_{idx}",
                        "points": [[point["t"], point[f"species_{idx}"]] for point in self._history],
                    }
                    for idx in range(1, self.species_count + 1)
                ],
            },
        }

    def _metrics_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Total abundance and diversity diagnostics for the Leibovich2022 community model.",
            "data": {
                "title": "Community Metrics",
                "series": [
                    {"name": "Total abundance", "points": [[point["t"], point["total_abundance"]] for point in self._history]},
                    {"name": "Species richness", "points": [[point["t"], point["richness"]] for point in self._history]},
                    {"name": "Shannon diversity", "points": [[point["t"], point["shannon_diversity"]] for point in self._history]},
                ],
            },
        }

    def _summary_visual(self) -> "VisualSpec":
        latest = self._history[-1]
        peak_total = max(point["total_abundance"] for point in self._history)
        min_total = min(point["total_abundance"] for point in self._history)
        peak_richness = max(point["richness"] for point in self._history)
        return {
            "render": "table",
            "description": "Competition, immigration, and diversity summary for the Leibovich2022 package.",
            "data": {
                "title": "Leibovich2022 Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["Species count", str(self.species_count)],
                    ["Competition overlap", f"{self.competition_overlap:.6g}"],
                    ["Immigration rate", f"{self.immigration_rate:.6g}"],
                    ["Final total abundance", f"{latest['total_abundance']:.6g}"],
                    ["Peak total abundance", f"{peak_total:.6g}"],
                    ["Minimum total abundance", f"{min_total:.6g}"],
                    ["Peak richness", f"{peak_richness:.6g}"],
                    ["Final Shannon diversity", f"{latest['shannon_diversity']:.6g}"],
                    ["Final dominant species", f"species_{int(latest['dominant_species_index'])}"],
                ],
            },
        }


Leibovich2022MultispeciesEcoCompetitionDescrModel2212080001Model = Leibovich2022CommunityModel

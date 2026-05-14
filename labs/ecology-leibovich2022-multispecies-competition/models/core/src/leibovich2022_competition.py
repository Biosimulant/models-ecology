# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Faithful stochastic multispecies competition model from Leibovich et al. (2022).

The upstream BioModels record for MODEL2212080001 distributes Python Gillespie
code, not SBML.  This Biosimulant package ports the upstream ``MultiLV``
propensity rules (gillespie_models.py) directly:

    birth propensity per species:  immi_rate + n_i * birth_rate
    death propensity per species:  n_i * (death_rate
        + (birth_rate - death_rate)
          * ((1 - comp_overlap) * n_i + comp_overlap * sum(n))
          / carry_capacity)

Immigration is uniform across species (no rarer-species weighting).
Competition coefficient is ``birth_rate - death_rate`` (not a separate
parameter), exactly matching the upstream ``MultiLV.propensity`` function.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, TYPE_CHECKING

import numpy as np
from biosim import StatefulBioModule
from biosim.signals import BioSignal, SignalSpec, coerce_float, scalar_or_record_input

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


class Leibovich2022CommunityModel(StatefulBioModule):
    """Faithful competition-immigration model matching upstream MultiLV propensities."""

    def __init__(
        self,
        species_count: int = 6,
        carrying_capacity: float = 100.0,
        birth_rate: float = 2.0,
        death_rate: float = 1.0,
        competition_overlap: float = 0.2,
        immigration_rate: float = 0.1,
        initial_abundance: float = 50.0,
        rng_seed: int = 7,
        integration_step: float = 0.05,
    ) -> None:
        if integration_step <= 0:
            raise ValueError("integration_step must be positive")
        if species_count < 2:
            raise ValueError("species_count must be at least 2")

        super().__init__(integration_step=integration_step, publish_on_zero_window=False)
        self.integration_step = float(integration_step)
        self.species_count = int(species_count)
        self.carrying_capacity = float(carrying_capacity)
        self.birth_rate = float(birth_rate)
        self.death_rate = float(death_rate)
        self.competition_overlap = float(competition_overlap)
        self.immigration_rate = float(immigration_rate)
        self.initial_abundance = float(initial_abundance)
        self.rng_seed = int(rng_seed)

        self._rng = np.random.default_rng(self.rng_seed)
        self._abundances = self._initial_state()

    def reset_state(self) -> None:
        self._rng = np.random.default_rng(self.rng_seed)
        self._abundances = self._initial_state()

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "carrying_capacity": scalar_or_record_input("individuals", "Habitat carrying capacity."),
            "birth_rate": scalar_or_record_input("1/time", "Per-capita birth rate."),
            "death_rate": scalar_or_record_input("1/time", "Per-capita baseline death rate."),
            "competition_overlap": scalar_or_record_input(
                "dimensionless", "Fraction of competition that is interspecific."
            ),
            "immigration_rate": scalar_or_record_input(
                "individuals/time", "Constant immigration rate per species."
            ),
            "initial_abundance": scalar_or_record_input("individuals", "Starting abundance per species."),
        }

    def _input_number(self, name: str) -> float | None:
        signal = self._input_overrides.get(name)
        if signal is None:
            return None
        return coerce_float(signal)

    def apply_overrides(self, *, reset_initial_state: bool) -> None:
        for attr in ("carrying_capacity", "birth_rate", "death_rate",
                     "competition_overlap", "immigration_rate"):
            value = self._input_number(attr)
            if value is not None and value >= 0.0:
                setattr(self, attr, value)

        value = self._input_number("initial_abundance")
        if value is not None and value >= 0.0:
            self.initial_abundance = value
            if reset_initial_state:
                self._abundances = self._initial_state()

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            'community_state': SignalSpec.record(schema={'species_abundances': 'json', 'total_abundance': 'json'}, emitted_unit='individuals', description='Species-resolved abundances for the stochastic competition community.'),
            'diversity_metrics': SignalSpec.record(schema={'richness': 'json', 'shannon_diversity': 'json', 'evenness': 'json', 'dominant_species': 'json'}, emitted_unit='dimensionless', description='Diversity and dominance metrics for the multispecies community.'),
            'visualisation_payload': SignalSpec.record(schema={'payload': 'json'}, description='Internal history payload for the sibling visualisation model.'),
        }

    def visualize(self) -> Optional[List["VisualSpec"]]:
        return None

    def _initial_state(self) -> np.ndarray:
        offsets = np.linspace(1.3, 0.7, self.species_count)
        return np.maximum(0, np.rint(self.initial_abundance * offsets)).astype(int)

    def step(self, dt: float) -> None:
        """Tau-leaping step with upstream MultiLV propensity formulas."""
        total_abundance = float(np.sum(self._abundances))
        next_state = self._abundances.copy()
        comp_coeff = self.birth_rate - self.death_rate

        for idx in range(self.species_count):
            n_i = float(self._abundances[idx])

            # Upstream birth propensity: immi_rate + n_i * birth_rate
            birth_prop = self.immigration_rate + n_i * self.birth_rate

            # Upstream death propensity: n_i * (death_rate + (b-d) * crowding / K)
            crowding = (
                (1.0 - self.competition_overlap) * n_i
                + self.competition_overlap * total_abundance
            )
            death_prop = n_i * (
                self.death_rate + comp_coeff * crowding / self.carrying_capacity
            )

            births = int(self._rng.poisson(max(0.0, birth_prop * dt)))
            deaths = int(self._rng.poisson(max(0.0, death_prop * dt)))
            next_state[idx] = max(0, int(n_i) + births - deaths)

        self._abundances = next_state

    def record_state(self, t: float) -> None:
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

    def output_payload(self, t: float) -> dict[str, Any]:
        latest = self._history[-1]
        species_abundances = {
            f"species_{idx}": int(self._abundances[idx - 1])
            for idx in range(1, self.species_count + 1)
        }
        return {
            "community_state": {
                    "species_abundances": species_abundances,
                    "total_abundance": int(latest["total_abundance"]),
                },
            "diversity_metrics": {
                    "richness": int(latest["richness"]),
                    "shannon_diversity": float(latest["shannon_diversity"]),
                    "evenness": float(latest["evenness"]),
                    "dominant_species": f"species_{int(latest['dominant_species_index'])}",
                },
            "visualisation_payload": {"payload": self._visualisation_payload()},
        }

    def _visualisation_payload(self) -> Dict[str, Any]:
        return {"species_count": self.species_count, "history": list(self._history)}

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

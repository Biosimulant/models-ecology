# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Reduced eco-genetic gene-drive suppression model grounded in Geci et al. (2022).

The upstream BioModels asset for MODEL2301120001 is a Julia implementation rather
than SBML. This package therefore exposes a curated reduced deterministic model
that keeps the ecological observables the paper is used for in Biosimulant:
population suppression, sex-ratio distortion, drive spread, and resistance.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

import biosim
from biosim.signals import BioSignal, SignalMetadata

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


class Geci2022GeneDriveModel(biosim.BioModule):
    """Reduced deterministic eco-genetic suppression model."""

    def __init__(
        self,
        initial_population: float = 12000.0,
        carrying_capacity: float = 15000.0,
        intrinsic_growth_rate: float = 0.85,
        base_mortality_rate: float = 0.22,
        initial_drive_frequency: float = 0.02,
        initial_resistance_frequency: float = 0.0,
        initial_male_fraction: float = 0.5,
        gene_drive_conversion_rate: float = 0.55,
        suppression_strength: float = 0.95,
        resistance_emergence_rate: float = 0.04,
        resistance_cost: float = 0.06,
        sex_ratio_bias: float = 0.22,
        min_dt: float = 0.1,
    ) -> None:
        if min_dt <= 0:
            raise ValueError("min_dt must be positive")
        if carrying_capacity <= 0 or initial_population < 0:
            raise ValueError("population sizes must be non-negative and capacity must be positive")

        self.min_dt = float(min_dt)
        self.initial_population = float(initial_population)
        self.carrying_capacity = float(carrying_capacity)
        self.intrinsic_growth_rate = float(intrinsic_growth_rate)
        self.base_mortality_rate = float(base_mortality_rate)
        self.initial_drive_frequency = float(initial_drive_frequency)
        self.initial_resistance_frequency = float(initial_resistance_frequency)
        self.initial_male_fraction = float(initial_male_fraction)
        self.gene_drive_conversion_rate = float(gene_drive_conversion_rate)
        self.suppression_strength = float(suppression_strength)
        self.resistance_emergence_rate = float(resistance_emergence_rate)
        self.resistance_cost = float(resistance_cost)
        self.sex_ratio_bias = float(sex_ratio_bias)

        self._time = 0.0
        self._population = self.initial_population
        self._drive_frequency = _clamp(self.initial_drive_frequency, 0.0, 1.0)
        self._resistance_frequency = _clamp(self.initial_resistance_frequency, 0.0, 1.0 - self._drive_frequency)
        self._male_fraction = _clamp(self.initial_male_fraction, 0.0, 1.0)
        self._history: List[Dict[str, float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return set()

    def outputs(self) -> Set[str]:
        return {"population_state", "gene_drive_metrics"}

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.reset()

    def reset(self) -> None:
        self._time = 0.0
        self._population = self.initial_population
        self._drive_frequency = _clamp(self.initial_drive_frequency, 0.0, 1.0)
        self._resistance_frequency = _clamp(self.initial_resistance_frequency, 0.0, 1.0 - self._drive_frequency)
        self._male_fraction = _clamp(self.initial_male_fraction, 0.0, 1.0)
        self._history = []
        self._outputs = {}

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
            self._population_visual(),
            self._genetics_visual(),
            self._summary_visual(),
        ]

    def _step(self, dt: float) -> None:
        wild_type = max(0.0, 1.0 - self._drive_frequency - self._resistance_frequency)
        effective_drive = self._drive_frequency * max(0.0, 1.0 - self._resistance_frequency)

        logistic_growth = self.intrinsic_growth_rate * (1.0 - self._population / self.carrying_capacity)
        suppression_load = self.suppression_strength * effective_drive
        net_growth = logistic_growth - self.base_mortality_rate - 0.6 * suppression_load
        self._population = max(0.0, self._population + dt * self._population * net_growth)

        drive_gain = self.gene_drive_conversion_rate * self._drive_frequency * wild_type
        drive_loss = self.resistance_emergence_rate * self._drive_frequency
        frequency_drag = 0.12 * self._drive_frequency * suppression_load
        self._drive_frequency = max(
            0.0,
            self._drive_frequency + dt * (drive_gain - drive_loss - frequency_drag),
        )

        resistance_gain = self.resistance_emergence_rate * self._drive_frequency * max(0.0, 1.0 - self._resistance_frequency)
        resistance_loss = self.resistance_cost * self._resistance_frequency * max(0.0, 1.0 - self._drive_frequency)
        self._resistance_frequency = max(0.0, self._resistance_frequency + dt * (resistance_gain - resistance_loss))

        combined = self._drive_frequency + self._resistance_frequency
        if combined > 1.0:
            scale = 1.0 / combined
            self._drive_frequency *= scale
            self._resistance_frequency *= scale

        self._male_fraction = _clamp(
            self.initial_male_fraction + self.sex_ratio_bias * effective_drive,
            0.5,
            0.98,
        )

    def _record_state(self, t: float) -> None:
        females = self._population * (1.0 - self._male_fraction)
        males = self._population * self._male_fraction
        suppression_ratio = 1.0 - (self._population / self.initial_population if self.initial_population > 0 else 0.0)
        self._history.append(
            {
                "t": float(t),
                "population": float(self._population),
                "adult_females": float(females),
                "adult_males": float(males),
                "drive_frequency": float(self._drive_frequency),
                "resistance_frequency": float(self._resistance_frequency),
                "male_fraction": float(self._male_fraction),
                "suppression_ratio": float(suppression_ratio),
            }
        )

    def _publish_outputs(self, t: float) -> None:
        latest = self._history[-1] if self._history else {
            "population": self._population,
            "adult_females": self._population * (1.0 - self._male_fraction),
            "adult_males": self._population * self._male_fraction,
            "drive_frequency": self._drive_frequency,
            "resistance_frequency": self._resistance_frequency,
            "male_fraction": self._male_fraction,
            "suppression_ratio": 1.0 - (self._population / self.initial_population if self.initial_population > 0 else 0.0),
        }
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "population_state": BioSignal(
                source=source_name,
                name="population_state",
                value={
                    "total_adults": float(latest["population"]),
                    "adult_females": float(latest["adult_females"]),
                    "adult_males": float(latest["adult_males"]),
                },
                time=float(t),
                metadata=SignalMetadata(
                    units="individuals",
                    description="Adult mosquito population partitioned into females and males.",
                    kind="state",
                ),
            ),
            "gene_drive_metrics": BioSignal(
                source=source_name,
                name="gene_drive_metrics",
                value={
                    "drive_frequency": float(latest["drive_frequency"]),
                    "resistance_frequency": float(latest["resistance_frequency"]),
                    "male_fraction": float(latest["male_fraction"]),
                    "suppression_ratio": float(latest["suppression_ratio"]),
                },
                time=float(t),
                metadata=SignalMetadata(
                    units=None,
                    description="Reduced eco-genetic metrics for drive spread, resistance, and suppression.",
                    kind="summary",
                ),
            ),
        }

    def _population_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Adult mosquito abundance through time, including sex-ratio distortion from the drive system.",
            "data": {
                "title": "Adult Population Suppression",
                "series": [
                    {"name": "Total adults", "points": [[p["t"], p["population"]] for p in self._history]},
                    {"name": "Adult females", "points": [[p["t"], p["adult_females"]] for p in self._history]},
                    {"name": "Adult males", "points": [[p["t"], p["adult_males"]] for p in self._history]},
                ],
            },
        }

    def _genetics_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Drive spread, resistance emergence, and male bias for the reduced Geci2022 eco-genetic model.",
            "data": {
                "title": "Gene-Drive Metrics",
                "series": [
                    {"name": "Drive frequency", "points": [[p["t"], p["drive_frequency"]] for p in self._history]},
                    {"name": "Resistance frequency", "points": [[p["t"], p["resistance_frequency"]] for p in self._history]},
                    {"name": "Male fraction", "points": [[p["t"], p["male_fraction"]] for p in self._history]},
                    {"name": "Suppression ratio", "points": [[p["t"], p["suppression_ratio"]] for p in self._history]},
                ],
            },
        }

    def _summary_visual(self) -> "VisualSpec":
        latest = self._history[-1]
        peak_drive = max(point["drive_frequency"] for point in self._history)
        peak_resistance = max(point["resistance_frequency"] for point in self._history)
        min_population = min(point["population"] for point in self._history)
        return {
            "render": "table",
            "description": "Final and extremal eco-genetic metrics for the reduced Geci2022 model.",
            "data": {
                "title": "Geci2022 Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["Initial adults", f"{self.initial_population:.6g}"],
                    ["Carrying capacity", f"{self.carrying_capacity:.6g}"],
                    ["Final adults", f"{latest['population']:.6g}"],
                    ["Minimum adults", f"{min_population:.6g}"],
                    ["Final female adults", f"{latest['adult_females']:.6g}"],
                    ["Peak drive frequency", f"{peak_drive:.6g}"],
                    ["Peak resistance frequency", f"{peak_resistance:.6g}"],
                    ["Final male fraction", f"{latest['male_fraction']:.6g}"],
                    ["Final suppression ratio", f"{latest['suppression_ratio']:.6g}"],
                ],
            },
        }


Geci2022Model2301120001Model = Geci2022GeneDriveModel

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

import biosim
import numpy as np
from biosim.signals import (AcceptedSignalProfile, ArraySignal, BioSignal, EventSignal, RecordSignal, ScalarSignal, SignalSpec)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


def _schema_type(value):
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int) and not isinstance(value, bool):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    return "json"


def _signal_value(signal):
    value = signal.value
    if isinstance(value, dict) and set(value.keys()) == {"payload"}:
        return value["payload"]
    return value


def _generic_input_spec(description=None):
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


def _make_signal(*, source, name, value, emitted_at, spec=None):
    if spec is None:
        if isinstance(value, dict):
            spec = SignalSpec.record(schema={str(key): _schema_type(item) for key, item in value.items()})
        elif isinstance(value, (list, tuple)):
            spec = SignalSpec.record(schema={"payload": "json"})
        else:
            spec = SignalSpec.scalar(dtype=_schema_type(value))

    if spec.signal_type == "scalar":
        return ScalarSignal(source=source, name=name, value=value, emitted_at=emitted_at, spec=spec)
    if spec.signal_type == "array":
        return ArraySignal(source=source, name=name, value=value, emitted_at=emitted_at, spec=spec)
    if spec.signal_type == "event":
        event_value = value
        if spec.schema is not None and not (isinstance(value, dict) and set(value.keys()) == set(spec.schema.keys())):
            event_value = {"payload": value}
        return EventSignal(source=source, name=name, value=event_value, emitted_at=emitted_at, spec=spec)

    record_value = value
    if not isinstance(value, dict) or set(value.keys()) != set((spec.schema or {}).keys()):
        record_value = {"payload": value}
    return RecordSignal(source=source, name=name, value=record_value, emitted_at=emitted_at, spec=spec)

class Leibovich2022CommunityModel(biosim.BioModule):
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

        self.integration_step = float(integration_step)
        self.species_count = int(species_count)
        self.carrying_capacity = float(carrying_capacity)
        self.birth_rate = float(birth_rate)
        self.death_rate = float(death_rate)
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

    def inputs(self) -> dict[str, SignalSpec]:
        return {}

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            'community_state': SignalSpec.record(schema={'species_abundances': 'json', 'total_abundance': 'json'}, emitted_unit='individuals', description='Species-resolved abundances for the stochastic competition community.'),
            'diversity_metrics': SignalSpec.record(schema={'richness': 'json', 'shannon_diversity': 'json', 'evenness': 'json', 'dominant_species': 'json'}, description='Diversity and dominance metrics for the multispecies community.'),
        }

    def advance_window(self, start: float, end: float) -> None:
        t = float(end)
        if t <= self._time:
            return

        current = self._time
        while current < t - 1e-12:
            dt = min(self.integration_step, t - current)
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
            "community_state": _make_signal(source=source_name, name="community_state", value={
                    "species_abundances": species_abundances,
                    "total_abundance": int(latest["total_abundance"]),
                }, emitted_at=float(t), spec=self.outputs().get("community_state")),
            "diversity_metrics": _make_signal(source=source_name, name="diversity_metrics", value={
                    "richness": int(latest["richness"]),
                    "shannon_diversity": float(latest["shannon_diversity"]),
                    "evenness": float(latest["evenness"]),
                    "dominant_species": f"species_{int(latest['dominant_species_index'])}",
                }, emitted_at=float(t), spec=self.outputs().get("diversity_metrics")),
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



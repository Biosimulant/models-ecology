# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Curated SBML packaging for Pfeiffer et al. (2001)."""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path
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

class Pfeiffer2001AtpProducingpathwaysCooperationcompetition(biosim.BioModule):
    """Expose pathway and resource observables from the Pfeiffer2001 SBML model."""

    _OBSERVABLES = {
        "S": "substrate_resource",
        "N1": "high_yield_population",
        "N2": "low_yield_population",
    }

    def __init__(self, model_path: str = "data/BIOMD0000000337.xml", integration_step: float = 0.1) -> None:
        if integration_step <= 0:
            raise ValueError("integration_step must be positive")
        self.integration_step = float(integration_step)
        self._model_path = Path(__file__).parent.parent / model_path
        self._time = 0.0
        self._runner: Any = None
        self._species_ids: List[str] = []
        self._species_units: Dict[str, Optional[str]] = {}
        self._history: List[Dict[str, float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self._runner = self._build_runner()
        self._species_ids = [sid for sid in self._runner.getFloatingSpeciesIds() if sid in self._OBSERVABLES]
        self._species_units = self._load_species_units()
        self._time = 0.0
        self._history = []
        self._outputs = {}

    def reset(self) -> None:
        if self._runner is not None and hasattr(self._runner, "reset"):
            self._runner.reset()
        self._time = 0.0
        self._history = []
        self._outputs = {}

    def inputs(self) -> dict[str, SignalSpec]:
        return {}

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            'resource_state': SignalSpec.record(schema={'substrate_resource': 'json'}, description='Shared substrate resource concentration from the Pfeiffer2001 SBML model.'),
            'community_state': SignalSpec.record(schema={'high_yield_population': 'json', 'low_yield_population': 'json', 'total_population': 'json'}, description='Population sizes for the two ATP-pathway strategies encoded by the model.'),
            'cooperation_metrics': SignalSpec.record(schema={'high_yield_fraction': 'json', 'low_yield_fraction': 'json', 'resource_per_biomass': 'json'}, description='Strategy fractions and resource intensity for the ATP-pathway competition model.'),
        }

    def advance_window(self, start: float, end: float) -> None:
        t = float(end)
        if self._runner is None:
            self.setup()
        if t <= self._time:
            return

        records = self._simulate_window(self._time, t)
        if records:
            self._history.extend(records)
            self._time = float(records[-1]["t"])
        else:
            self._time = float(t)

        self._publish_outputs(self._time)

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[List["VisualSpec"]]:
        if not self._history:
            return None
        return [
            self._trajectory_visual(),
            self._fraction_visual(),
            self._summary_visual(),
        ]

    def _build_runner(self) -> Any:
        import tellurium as te

        return te.loadSBMLModel(str(self._model_path))

    def _load_species_units(self) -> Dict[str, Optional[str]]:
        tree = ET.parse(self._model_path)
        root = tree.getroot()
        ns = {"sbml": root.tag.split("}")[0].strip("{")}
        species_units: Dict[str, Optional[str]] = {}
        for species in root.findall(".//sbml:species", ns):
            species_units[species.attrib["id"]] = species.attrib.get("substanceUnits") or species.attrib.get("units")
        return species_units

    def _simulate_window(self, start: float, end: float) -> List[Dict[str, float]]:
        step_count = max(2, int(math.ceil((end - start) / self.integration_step)) + 1)
        selections = ["time", *self._species_ids]
        result = self._runner.simulate(start, end, step_count, selections=selections)
        matrix = np.asarray(result, dtype=float)
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)

        records: List[Dict[str, float]] = []
        for row_index, row in enumerate(matrix):
            if row_index == 0 and self._history and abs(float(row[0]) - self._history[-1]["t"]) < 1e-12:
                continue
            records.append(
                {
                    "t": float(row[0]),
                    "S": float(row[1]),
                    "N1": float(row[2]),
                    "N2": float(row[3]),
                }
            )
        return records

    def _publish_outputs(self, t: float) -> None:
        latest = self._history[-1]
        total_population = latest["N1"] + latest["N2"]
        high_yield_fraction = latest["N1"] / total_population if total_population > 0 else 0.0
        low_yield_fraction = latest["N2"] / total_population if total_population > 0 else 0.0
        resource_per_biomass = latest["S"] / total_population if total_population > 0 else 0.0
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "resource_state": _make_signal(source=source_name, name="resource_state", value={"substrate_resource": latest["S"]}, emitted_at=float(t), spec=self.outputs().get("resource_state") if 'self' in locals() else None),
            "community_state": _make_signal(source=source_name, name="community_state", value={
                    "high_yield_population": latest["N1"],
                    "low_yield_population": latest["N2"],
                    "total_population": total_population,
                }, emitted_at=float(t), spec=self.outputs().get("community_state") if 'self' in locals() else None),
            "cooperation_metrics": _make_signal(source=source_name, name="cooperation_metrics", value={
                    "high_yield_fraction": high_yield_fraction,
                    "low_yield_fraction": low_yield_fraction,
                    "resource_per_biomass": resource_per_biomass,
                }, emitted_at=float(t), spec=self.outputs().get("cooperation_metrics") if 'self' in locals() else None),
        }

    def _trajectory_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Substrate and strategy-specific biomass trajectories from the Pfeiffer2001 ATP-pathway model.",
            "data": {
                "title": "Resource and Population Trajectories",
                "series": [
                    {"name": "Substrate resource", "points": [[point["t"], point["S"]] for point in self._history]},
                    {"name": "High-yield population", "points": [[point["t"], point["N1"]] for point in self._history]},
                    {"name": "Low-yield population", "points": [[point["t"], point["N2"]] for point in self._history]},
                ],
            },
        }

    def _fraction_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Population-fraction view for the two ATP-producing strategies, plus substrate per unit biomass.",
            "data": {
                "title": "Strategy Fractions",
                "series": [
                    {
                        "name": "High-yield fraction",
                        "points": [
                            [point["t"], point["N1"] / (point["N1"] + point["N2"]) if (point["N1"] + point["N2"]) > 0 else 0.0]
                            for point in self._history
                        ],
                    },
                    {
                        "name": "Low-yield fraction",
                        "points": [
                            [point["t"], point["N2"] / (point["N1"] + point["N2"]) if (point["N1"] + point["N2"]) > 0 else 0.0]
                            for point in self._history
                        ],
                    },
                    {
                        "name": "Resource per biomass",
                        "points": [
                            [point["t"], point["S"] / (point["N1"] + point["N2"]) if (point["N1"] + point["N2"]) > 0 else 0.0]
                            for point in self._history
                        ],
                    },
                ],
            },
        }

    def _summary_visual(self) -> "VisualSpec":
        latest = self._history[-1]
        total_population = latest["N1"] + latest["N2"]
        peak_total = max(point["N1"] + point["N2"] for point in self._history)
        min_resource = min(point["S"] for point in self._history)
        return {
            "render": "table",
            "description": "Summary statistics for substrate depletion and pathway competition in the Pfeiffer2001 model.",
            "data": {
                "title": "Pfeiffer2001 Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["Final substrate", f"{latest['S']:.6g}"],
                    ["Minimum substrate", f"{min_resource:.6g}"],
                    ["Final high-yield population", f"{latest['N1']:.6g}"],
                    ["Final low-yield population", f"{latest['N2']:.6g}"],
                    ["Final total population", f"{total_population:.6g}"],
                    ["Peak total population", f"{peak_total:.6g}"],
                    ["Final high-yield fraction", f"{(latest['N1'] / total_population) if total_population > 0 else 0.0:.6g}"],
                    ["Final low-yield fraction", f"{(latest['N2'] / total_population) if total_population > 0 else 0.0:.6g}"],
                ],
            },
        }


Pfeiffer2001AtpProducingpathwaysCooperationcoBiomd0000000337Model = Pfeiffer2001AtpProducingpathwaysCooperationcompetition

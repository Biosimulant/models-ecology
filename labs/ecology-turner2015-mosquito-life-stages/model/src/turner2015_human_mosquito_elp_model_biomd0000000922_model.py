# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Curated SBML packaging for the Turner2015 mosquito life-stage model."""
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

class SbmlTurner2015HumanMosquitoElpModel(biosim.BioModule):
    """Expose life-stage ecology observables from the Turner2015 SBML model."""

    _STAGES = {
        "Population_of_Eggs": "eggs",
        "Population_of_Larvae": "larvae",
        "Population_of_Pupae": "pupae",
    }

    def __init__(self, model_path: str = "data/BIOMD0000000922.xml", integration_step: float = 0.1) -> None:
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
        self._species_ids = [sid for sid in self._runner.getFloatingSpeciesIds() if sid in self._STAGES]
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
            'life_stage_state': SignalSpec.record(schema={'eggs': 'json', 'larvae': 'json', 'pupae': 'json', 'total_immature_population': 'json'}, description='Immature mosquito life-stage abundances from the Turner2015 model.'),
            'population_metrics': SignalSpec.record(schema={'egg_fraction': 'json', 'larval_fraction': 'json', 'pupal_fraction': 'json', 'total_immature_population': 'json'}, description='Stage fractions and total immature abundance for the Turner2015 mosquito life-stage model.'),
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
            self._life_stage_visual(),
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
                    "Population_of_Eggs": float(row[1]),
                    "Population_of_Larvae": float(row[2]),
                    "Population_of_Pupae": float(row[3]),
                }
            )
        return records

    def _publish_outputs(self, t: float) -> None:
        latest = self._history[-1]
        eggs = latest["Population_of_Eggs"]
        larvae = latest["Population_of_Larvae"]
        pupae = latest["Population_of_Pupae"]
        total_immature = eggs + larvae + pupae
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "life_stage_state": _make_signal(source=source_name, name="life_stage_state", value={
                    "eggs": eggs,
                    "larvae": larvae,
                    "pupae": pupae,
                    "total_immature_population": total_immature,
                }, emitted_at=float(t), spec=self.outputs().get("life_stage_state") if 'self' in locals() else None),
            "population_metrics": _make_signal(source=source_name, name="population_metrics", value={
                    "egg_fraction": eggs / total_immature if total_immature > 0 else 0.0,
                    "larval_fraction": larvae / total_immature if total_immature > 0 else 0.0,
                    "pupal_fraction": pupae / total_immature if total_immature > 0 else 0.0,
                    "total_immature_population": total_immature,
                }, emitted_at=float(t), spec=self.outputs().get("population_metrics") if 'self' in locals() else None),
        }

    def _life_stage_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Egg, larval, and pupal mosquito abundances through time.",
            "data": {
                "title": "Mosquito Immature Life Stages",
                "series": [
                    {
                        "name": "Eggs",
                        "points": [[point["t"], point["Population_of_Eggs"]] for point in self._history],
                    },
                    {
                        "name": "Larvae",
                        "points": [[point["t"], point["Population_of_Larvae"]] for point in self._history],
                    },
                    {
                        "name": "Pupae",
                        "points": [[point["t"], point["Population_of_Pupae"]] for point in self._history],
                    },
                ],
            },
        }

    def _fraction_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Stage-fraction diagnostics for the Turner2015 immature mosquito population.",
            "data": {
                "title": "Life-Stage Fractions",
                "series": [
                    {
                        "name": "Egg fraction",
                        "points": [
                            [
                                point["t"],
                                point["Population_of_Eggs"]
                                / (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"])
                                if (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"]) > 0
                                else 0.0,
                            ]
                            for point in self._history
                        ],
                    },
                    {
                        "name": "Larval fraction",
                        "points": [
                            [
                                point["t"],
                                point["Population_of_Larvae"]
                                / (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"])
                                if (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"]) > 0
                                else 0.0,
                            ]
                            for point in self._history
                        ],
                    },
                    {
                        "name": "Pupal fraction",
                        "points": [
                            [
                                point["t"],
                                point["Population_of_Pupae"]
                                / (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"])
                                if (point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"]) > 0
                                else 0.0,
                            ]
                            for point in self._history
                        ],
                    },
                ],
            },
        }

    def _summary_visual(self) -> "VisualSpec":
        latest = self._history[-1]
        totals = [
            point["Population_of_Eggs"] + point["Population_of_Larvae"] + point["Population_of_Pupae"]
            for point in self._history
        ]
        return {
            "render": "table",
            "description": "Summary statistics for the mosquito life-stage composition encoded by Turner2015.",
            "data": {
                "title": "Turner2015 Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["Final eggs", f"{latest['Population_of_Eggs']:.6g}"],
                    ["Final larvae", f"{latest['Population_of_Larvae']:.6g}"],
                    ["Final pupae", f"{latest['Population_of_Pupae']:.6g}"],
                    ["Final total immature", f"{totals[-1]:.6g}"],
                    ["Peak total immature", f"{max(totals):.6g}"],
                    ["Minimum total immature", f"{min(totals):.6g}"],
                ],
            },
        }


Turner2015HumanMosquitoElpModelBiomd0000000922Model = SbmlTurner2015HumanMosquitoElpModel

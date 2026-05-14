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
from biosim.signals import unwrap_payload as _signal_value
from biosim.signals import make_signal as _make_signal

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



def _generic_input_spec(description=None):
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


class Turner2015MosquitoLifeStagesModel(biosim.BioModule):
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
        self._input_overrides: Dict[str, BioSignal] = {}
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

    @staticmethod
    def _scalar_input_spec(unit: str, description: str) -> SignalSpec:
        return SignalSpec.scalar(
            dtype="float64",
            accepted_profiles=(
                AcceptedSignalProfile(signal_type="scalar", dtype="float64",
                                     accepted_units=(unit,), description=description),
                AcceptedSignalProfile(signal_type="record", schema={"payload": "json"},
                                     description=description),
            ),
            description=description,
        )

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "integration_step": self._scalar_input_spec("day", "ODE solver integration step size."),
        }

    def set_inputs(self, inputs: dict[str, BioSignal]) -> None:
        self._input_overrides = dict(inputs or {})
        self._apply_input_overrides()

    def _input_number(self, name: str) -> float | None:
        signal = self._input_overrides.get(name)
        if signal is None:
            return None
        value = _signal_value(signal)
        if isinstance(value, dict):
            for key in ("value", "count", "payload"):
                if key in value:
                    value = value[key]
                    break
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _apply_input_overrides(self) -> None:
        value = self._input_number("integration_step")
        if value is not None and value > 0:
            self.integration_step = value

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            'life_stage_state': SignalSpec.record(schema={'eggs': 'json', 'larvae': 'json', 'pupae': 'json', 'total_immature_population': 'json'}, emitted_unit='individuals', description='Immature mosquito life-stage abundances from the Turner2015 model.'),
            'population_metrics': SignalSpec.record(schema={'egg_fraction': 'json', 'larval_fraction': 'json', 'pupal_fraction': 'json', 'total_immature_population': 'json'}, emitted_unit='fraction', description='Stage fractions and total immature abundance for the Turner2015 mosquito life-stage model.'),
            'visualisation_payload': SignalSpec.record(schema={'payload': 'json'}, description='Internal history payload for the sibling visualisation model.'),
        }

    def advance_window(self, start: float, end: float, inputs: dict[str, BioSignal] | None = None) -> None:
        if inputs:
            self.set_inputs(inputs)
        else:
            self._apply_input_overrides()

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
        return None

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
            "visualisation_payload": _make_signal(
                source=source_name,
                name="visualisation_payload",
                value={"payload": self._visualisation_payload()},
                emitted_at=float(t),
                spec=self.outputs().get("visualisation_payload"),
            ),
        }

    def _visualisation_payload(self) -> Dict[str, Any]:
        return {"history": list(self._history)}

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


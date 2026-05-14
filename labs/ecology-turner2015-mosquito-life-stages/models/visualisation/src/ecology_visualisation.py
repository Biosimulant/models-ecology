# SPDX-FileCopyrightText: 2026-present Biosimulant Team
# SPDX-License-Identifier: Apache-2.0
"""Dedicated visualisation model for ecology labs."""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional

from biosim import BioModule
from biosim.signals import AcceptedSignalProfile, BioSignal, SignalSpec
from biosim.signals import unwrap_payload as _signal_value



def _record_input_spec(description: str) -> SignalSpec:
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


class EcologyVisualisationModel(BioModule):
    def __init__(self, integration_step: float = 0.1, source_alias: str = "core", mode: str = "lotka_volterra", lab_title: str = "Ecology Lab") -> None:
        self.integration_step = float(integration_step)
        self.source_alias = source_alias
        self.mode = mode
        self.lab_title = lab_title
        self._inputs: Dict[str, BioSignal] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        return {f"{self.source_alias}_visualisation_payload": _record_input_spec("Internal payload from the sibling core model.")}

    def outputs(self) -> dict[str, SignalSpec]:
        return {}

    def reset(self) -> None:
        self._inputs = {}

    def set_inputs(self, signals: dict[str, BioSignal]) -> None:
        self._inputs.update(signals or {})

    def advance_window(self, start: float | None = None, end: float | None = None, inputs: dict[str, BioSignal] | None = None) -> dict[str, BioSignal]:
        if inputs:
            self.set_inputs(inputs)
        return {}

    def get_outputs(self) -> dict[str, BioSignal]:
        return {}

    def visualize(self) -> Optional[list[dict[str, Any]]]:
        payload = _signal_value(self._inputs.get(f"{self.source_alias}_visualisation_payload"))
        if isinstance(payload, dict) and set(payload.keys()) == {"payload"}:
            payload = payload["payload"]
        if not isinstance(payload, Mapping):
            return None
        history = payload.get("history")
        if not isinstance(history, list) or not history:
            return None
        if self.mode == "lotka_volterra":
            return self._lotka_visuals(payload, history)
        if self.mode == "pfeiffer":
            return self._pfeiffer_visuals(payload, history)
        if self.mode == "turner":
            return self._turner_visuals(payload, history)
        if self.mode == "leibovich":
            return self._leibovich_visuals(payload, history)
        if self.mode == "gene_drive":
            return self._gene_drive_visuals(payload, history)
        return self._rosenzweig_visuals(payload, history)

    def _lotka_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        params = payload.get("parameters", {}) if isinstance(payload.get("parameters"), Mapping) else {}
        prey_name = params.get("prey_name", "Prey")
        predator_name = params.get("predator_name", "Predator")
        prey_values = [float(point.get("prey", 0.0)) for point in history]
        predator_values = [float(point.get("predator", 0.0)) for point in history]
        alpha = float(params.get("alpha", 0.0))
        beta = float(params.get("beta", 0.0))
        gamma = float(params.get("gamma", 0.0))
        delta = float(params.get("delta", 0.0))
        rows = [
            ["prey growth rate alpha", f"{alpha:.6g}", "1/day"],
            ["predation rate beta", f"{beta:.6g}", "1/(count*day)"],
            ["predator mortality gamma", f"{gamma:.6g}", "1/day"],
            ["predator reproduction delta", f"{delta:.6g}", "1/(count*day)"],
            ["peak prey", f"{max(prey_values):.6g}", "count"],
            ["peak predator", f"{max(predator_values):.6g}", "count"],
            ["prey extinction time", str(payload.get("prey_extinction_time")), "day"],
            ["predator extinction time", str(payload.get("predator_extinction_time")), "day"],
        ]
        return [
            {"render": "timeseries", "description": "Predator-prey population trajectories.", "data": {"title": "Population Trajectories", "x_unit": "day", "y_unit": "count", "series": [{"name": prey_name, "points": [[float(p.get('t', 0.0)), float(p.get('prey', 0.0))] for p in history]}, {"name": predator_name, "points": [[float(p.get('t', 0.0)), float(p.get('predator', 0.0))] for p in history]}]}},
            {"render": "timeseries", "description": "Phase trajectory expressed as predator count against prey count.", "data": {"title": "Phase Trajectory", "x_unit": "count", "y_unit": "count", "series": [{"name": f"{prey_name} vs {predator_name}", "points": [[float(p.get('prey', 0.0)), float(p.get('predator', 0.0))] for p in history]}]}},
            {"render": "table", "description": "Lotka-Volterra summary diagnostics.", "data": {"title": "Lotka-Volterra Summary", "columns": ["Metric", "Value", "Unit"], "rows": rows}},
            {"render": "timeseries", "description": "Invariant and drift diagnostics across the run.", "data": {"title": "Invariant Audit", "x_unit": "day", "series": [{"name": "Invariant", "points": [[float(p.get('t', 0.0)), float(p.get('invariant', 0.0))] for p in history]}, {"name": "Drift from initial", "points": [[float(p.get('t', 0.0)), float(p.get('drift', 0.0))] for p in history]}]}},
        ]

    def _pfeiffer_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        latest = history[-1]
        totals = [float(p.get("N1", 0.0)) + float(p.get("N2", 0.0)) for p in history]
        return [
            {"render": "timeseries", "description": "Substrate and pathway-specific population trajectories.", "data": {"title": "Resource and Population Trajectories", "series": [{"name": "Substrate resource", "points": [[float(p.get('t', 0.0)), float(p.get('S', 0.0))] for p in history]}, {"name": "High-yield population", "points": [[float(p.get('t', 0.0)), float(p.get('N1', 0.0))] for p in history]}, {"name": "Low-yield population", "points": [[float(p.get('t', 0.0)), float(p.get('N2', 0.0))] for p in history]}]}},
            {"render": "timeseries", "description": "Strategy fractions and resource intensity through time.", "data": {"title": "Strategy Fractions", "series": [{"name": "High-yield fraction", "points": [[float(p.get('t', 0.0)), float(p.get('N1', 0.0)) / max(float(p.get('N1', 0.0)) + float(p.get('N2', 0.0)), 1e-12)] for p in history]}, {"name": "Low-yield fraction", "points": [[float(p.get('t', 0.0)), float(p.get('N2', 0.0)) / max(float(p.get('N1', 0.0)) + float(p.get('N2', 0.0)), 1e-12)] for p in history]}, {"name": "Resource per biomass", "points": [[float(p.get('t', 0.0)), float(p.get('S', 0.0)) / max(float(p.get('N1', 0.0)) + float(p.get('N2', 0.0)), 1e-12)] for p in history]}]}},
            {"render": "table", "description": "Summary statistics for ATP-pathway competition.", "data": {"title": "Pfeiffer2001 Summary", "columns": ["Metric", "Value"], "rows": [["Final substrate", f"{float(latest.get('S', 0.0)):.6g}"], ["Minimum substrate", f"{min(float(p.get('S', 0.0)) for p in history):.6g}"], ["Final total population", f"{totals[-1]:.6g}"], ["Peak total population", f"{max(totals):.6g}"]] }},
        ]

    def _turner_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        latest = history[-1]
        totals = [float(p.get("Population_of_Eggs", 0.0)) + float(p.get("Population_of_Larvae", 0.0)) + float(p.get("Population_of_Pupae", 0.0)) for p in history]
        return [
            {"render": "timeseries", "description": "Mosquito immature life stages over time.", "data": {"title": "Mosquito Immature Life Stages", "series": [{"name": "Eggs", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Eggs', 0.0))] for p in history]}, {"name": "Larvae", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Larvae', 0.0))] for p in history]}, {"name": "Pupae", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Pupae', 0.0))] for p in history]}]}},
            {"render": "timeseries", "description": "Life-stage fractions across the run.", "data": {"title": "Life-Stage Fractions", "series": [{"name": "Egg fraction", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Eggs', 0.0)) / max(float(p.get('Population_of_Eggs', 0.0)) + float(p.get('Population_of_Larvae', 0.0)) + float(p.get('Population_of_Pupae', 0.0)), 1e-12)] for p in history]}, {"name": "Larval fraction", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Larvae', 0.0)) / max(float(p.get('Population_of_Eggs', 0.0)) + float(p.get('Population_of_Larvae', 0.0)) + float(p.get('Population_of_Pupae', 0.0)), 1e-12)] for p in history]}, {"name": "Pupal fraction", "points": [[float(p.get('t', 0.0)), float(p.get('Population_of_Pupae', 0.0)) / max(float(p.get('Population_of_Eggs', 0.0)) + float(p.get('Population_of_Larvae', 0.0)) + float(p.get('Population_of_Pupae', 0.0)), 1e-12)] for p in history]}]}},
            {"render": "table", "description": "Summary statistics for the life-stage composition.", "data": {"title": "Turner2015 Summary", "columns": ["Metric", "Value"], "rows": [["Final eggs", f"{float(latest.get('Population_of_Eggs', 0.0)):.6g}"], ["Final larvae", f"{float(latest.get('Population_of_Larvae', 0.0)):.6g}"], ["Final pupae", f"{float(latest.get('Population_of_Pupae', 0.0)):.6g}"], ["Peak total immature", f"{max(totals):.6g}"]] }},
        ]

    def _leibovich_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        species_count = int(payload.get("species_count", 0))
        latest = history[-1]
        return [
            {"render": "timeseries", "description": "Species-resolved abundance trajectories.", "data": {"title": "Community Abundance Trajectories", "series": [{"name": f"species_{idx}", "points": [[float(p.get('t', 0.0)), float(p.get(f'species_{idx}', 0.0))] for p in history]} for idx in range(1, species_count + 1)]}},
            {"render": "timeseries", "description": "Community-level diversity diagnostics.", "data": {"title": "Community Metrics", "series": [{"name": "Total abundance", "points": [[float(p.get('t', 0.0)), float(p.get('total_abundance', 0.0))] for p in history]}, {"name": "Species richness", "points": [[float(p.get('t', 0.0)), float(p.get('richness', 0.0))] for p in history]}, {"name": "Shannon diversity", "points": [[float(p.get('t', 0.0)), float(p.get('shannon_diversity', 0.0))] for p in history]}]}},
            {"render": "table", "description": "Competition and diversity summary.", "data": {"title": "Leibovich2022 Summary", "columns": ["Metric", "Value"], "rows": [["Species count", str(species_count)], ["Final total abundance", f"{float(latest.get('total_abundance', 0.0)):.6g}"], ["Peak richness", f"{max(float(p.get('richness', 0.0)) for p in history):.6g}"], ["Dominant species", f"species_{int(float(latest.get('dominant_species_index', 1.0)))}"]]}},
        ]

    def _gene_drive_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        latest = history[-1]
        return [
            {"render": "timeseries", "description": "Adult mosquito abundance and sex split.", "data": {"title": "Adult Population Suppression", "series": [{"name": "Total adults", "points": [[float(p.get('t', 0.0)), float(p.get('total_adults', 0.0))] for p in history]}, {"name": "Adult females", "points": [[float(p.get('t', 0.0)), float(p.get('adult_females', 0.0))] for p in history]}, {"name": "Adult males", "points": [[float(p.get('t', 0.0)), float(p.get('adult_males', 0.0))] for p in history]}]}},
            {"render": "timeseries", "description": "Drive spread and resistance metrics.", "data": {"title": "Gene-Drive Metrics", "series": [{"name": "Drive frequency", "points": [[float(p.get('t', 0.0)), float(p.get('drive_frequency', 0.0))] for p in history]}, {"name": "Resistance frequency", "points": [[float(p.get('t', 0.0)), float(p.get('resistance_frequency', 0.0))] for p in history]}, {"name": "Male fraction", "points": [[float(p.get('t', 0.0)), float(p.get('male_fraction', 0.0))] for p in history]}, {"name": "Suppression ratio", "points": [[float(p.get('t', 0.0)), float(p.get('suppression_ratio', 0.0))] for p in history]}]}},
            {"render": "table", "description": "Eco-genetic summary metrics.", "data": {"title": "Geci2022 Summary", "columns": ["Metric", "Value"], "rows": [["Initial total", f"{float(payload.get('initial_total', 0.0)):.6g}"], ["Final total adults", f"{float(latest.get('total_adults', 0.0)):.6g}"], ["Peak drive frequency", f"{max(float(p.get('drive_frequency', 0.0)) for p in history):.6g}"], ["Final suppression ratio", f"{float(latest.get('suppression_ratio', 0.0)):.6g}"]] }},
        ]

    def _rosenzweig_visuals(self, payload: Mapping[str, Any], history: list[dict[str, Any]]) -> list[dict[str, Any]]:
        prey_label = str(payload.get("prey_label", "Prey"))
        predator_label = str(payload.get("predator_label", "Predator"))
        equilibrium = payload.get("equilibrium_summary", {}) if isinstance(payload.get("equilibrium_summary"), Mapping) else {}
        stability = payload.get("stability_summary", {}) if isinstance(payload.get("stability_summary"), Mapping) else {}
        risk = payload.get("extinction_risk", {}) if isinstance(payload.get("extinction_risk"), Mapping) else {}
        return [
            {"render": "timeseries", "description": "Population trajectories for the configured Rosenzweig-MacArthur scenario.", "data": {"title": "Population Trajectories", "x_unit": payload.get("time_unit", "day"), "y_unit": "count", "series": [{"name": prey_label, "points": [[float(p.get('t', 0.0)), float(p.get('prey', 0.0))] for p in history]}, {"name": predator_label, "points": [[float(p.get('t', 0.0)), float(p.get('predator', 0.0))] for p in history]}]}},
            {"render": "timeseries", "description": "Phase trajectory expressed as predator count against prey count.", "data": {"title": "Phase Trajectory", "x_unit": "count", "y_unit": "count", "series": [{"name": f"{prey_label} vs {predator_label}", "points": [[float(p.get('prey', 0.0)), float(p.get('predator', 0.0))] for p in history]}]}},
            {"render": "table", "description": "Equilibrium, stability, and risk summary.", "data": {"title": "Ecology State and Risk", "columns": ["Metric", "Value"], "rows": [["Prey equilibrium", str(equilibrium.get('prey_equilibrium'))], ["Predator equilibrium", str(equilibrium.get('predator_equilibrium'))], ["Regime", str(stability.get('regime'))], ["Prey extinction risk", str(risk.get('prey'))], ["Predator extinction risk", str(risk.get('predator'))]]}},
        ]

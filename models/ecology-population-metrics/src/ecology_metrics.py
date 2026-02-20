# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Ecology metrics: summary statistics for ecological simulations."""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim import BioWorld
    from biosim.visuals import VisualSpec

from biosim import BioModule
from biosim.signals import BioSignal, SignalMetadata


class EcologyMetrics(BioModule):
    """Compute summary statistics for ecological simulations.

    Produces a table with metrics like:
    - Total population across all species
    - Species diversity (Shannon index)
    - Population stability (coefficient of variation)
    - Extinction events
    """

    def __init__(self, min_dt: float = 1.0) -> None:
        self.min_dt = min_dt
        self._populations: Dict[str, List[int]] = {}
        self._extinctions: Dict[str, float] = {}  # species -> extinction time
        self._t_start: Optional[float] = None
        self._t_end: float = 0.0
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"population_state"}

    def outputs(self) -> Set[str]:
        return {"metrics"}

    def reset(self) -> None:
        self._populations = {}
        self._extinctions = {}
        self._t_start = None
        self._t_end = 0.0
        self._outputs = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        signal = signals.get("population_state")
        if signal is None or not isinstance(signal.value, dict):
            return
        species = str(signal.value.get("species", "Unknown"))
        count = int(signal.value.get("count", 0))
        t = float(signal.value.get("t", signal.time))

        if self._t_start is None:
            self._t_start = t
        self._t_end = t

        if species not in self._populations:
            self._populations[species] = []

        self._populations[species].append(count)

        # Track extinctions
        if count == 0 and species not in self._extinctions:
            self._extinctions[species] = t

    def advance_to(self, t: float) -> None:
        # Emit incremental metrics as a state signal for persistence.
        self._t_end = max(self._t_end, float(t))
        duration = self._t_end - (self._t_start or 0.0)

        final_total = sum((h[-1] if h else 0) for h in self._populations.values())
        n_species = len(self._populations)
        n_extinct = len(self._extinctions)
        n_surviving = n_species - n_extinct
        shannon = self._compute_shannon_diversity()

        cvs = []
        for history in self._populations.values():
            cv = self._compute_cv(history)
            if cv is not None:
                cvs.append(cv)
        avg_cv = sum(cvs) / len(cvs) if cvs else 0.0

        source = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "metrics": BioSignal(
                source=source,
                name="metrics",
                value={
                    "t": float(t),
                    "duration": float(duration),
                    "n_species": int(n_species),
                    "n_surviving": int(n_surviving),
                    "n_extinct": int(n_extinct),
                    "final_total_population": int(final_total),
                    "shannon_diversity": float(shannon),
                    "avg_population_cv": float(avg_cv),
                    "extinctions": {k: float(v) for k, v in self._extinctions.items()},
                },
                time=float(t),
                metadata=SignalMetadata(description="Ecology summary metrics", kind="state"),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def _compute_shannon_diversity(self) -> float:
        """Compute Shannon diversity index from final populations."""
        import math

        final_pops = {}
        for species, history in self._populations.items():
            if history:
                final_pops[species] = history[-1]

        total = sum(final_pops.values())
        if total <= 0:
            return 0.0

        H = 0.0
        for count in final_pops.values():
            if count > 0:
                p = count / total
                H -= p * math.log(p)

        return H

    def _compute_cv(self, values: List[int]) -> Optional[float]:
        """Compute coefficient of variation."""
        if len(values) < 2:
            return None

        mean = sum(values) / len(values)
        if mean <= 0:
            return None

        variance = sum((v - mean) ** 2 for v in values) / len(values)
        std = variance ** 0.5
        return std / mean

    def visualize(self) -> Optional["VisualSpec"]:
        """Generate ecology metrics table."""
        duration = self._t_end - (self._t_start or 0.0)

        # Compute metrics
        final_total = sum(
            history[-1] if history else 0
            for history in self._populations.values()
        )

        n_species = len(self._populations)
        n_extinct = len(self._extinctions)
        n_surviving = n_species - n_extinct

        shannon = self._compute_shannon_diversity()

        # Average CV across species (population stability)
        cvs = []
        for history in self._populations.values():
            cv = self._compute_cv(history)
            if cv is not None:
                cvs.append(cv)
        avg_cv = sum(cvs) / len(cvs) if cvs else 0.0

        # Peak total population
        max_totals = []
        if self._populations:
            first_species = list(self._populations.keys())[0]
            n_steps = len(self._populations[first_species])
            for i in range(n_steps):
                total = sum(
                    history[i] if i < len(history) else 0
                    for history in self._populations.values()
                )
                max_totals.append(total)
        peak_pop = max(max_totals) if max_totals else 0

        rows = [
            ["Duration (time units)", f"{duration:.2f}"],
            ["Number of Species", str(n_species)],
            ["Surviving Species", str(n_surviving)],
            ["Extinctions", str(n_extinct)],
            ["Final Total Population", str(final_total)],
            ["Peak Total Population", str(peak_pop)],
            ["Shannon Diversity Index", f"{shannon:.3f}"],
            ["Avg Population CV", f"{avg_cv:.3f}"],
        ]

        # Add extinction details
        if self._extinctions:
            for species, ext_time in sorted(self._extinctions.items()):
                rows.append([f"  {species} Extinct at", f"t={ext_time:.2f}"])

        return {
            "render": "table",
            "data": {
                "columns": ["Metric", "Value"],
                "rows": rows,
            },
        }

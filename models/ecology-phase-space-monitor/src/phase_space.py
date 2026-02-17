# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Phase space monitor: 2D phase space plot for two-species dynamics."""
from __future__ import annotations

import base64
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from bsim import BioWorld
    from bsim.visuals import VisualSpec

from bsim import BioModule
from bsim.signals import BioSignal, SignalMetadata


class PhaseSpaceMonitor(BioModule):
    """Creates a phase space plot of two populations.

    Useful for visualizing predator-prey dynamics in 2D space.

    Parameters:
        x_species: Name of species for X axis.
        y_species: Name of species for Y axis.
        max_points: Maximum points to store.
    """

    def __init__(
        self,
        x_species: str = "Prey",
        y_species: str = "Predator",
        max_points: int = 5000,
        min_dt: float = 1.0,
    ) -> None:
        self.min_dt = min_dt
        self.x_species = x_species
        self.y_species = y_species
        self.max_points = max_points

        self._x_values: List[int] = []
        self._y_values: List[int] = []
        self._current_x: int = 0
        self._current_y: int = 0
        self._time: float = 0.0
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"population_state"}

    def outputs(self) -> Set[str]:
        return {"phase_point"}

    def reset(self) -> None:
        self._x_values = []
        self._y_values = []
        self._current_x = 0
        self._current_y = 0
        self._time = 0.0
        self._outputs = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        signal = signals.get("population_state")
        if signal is None or not isinstance(signal.value, dict):
            return
        species = str(signal.value.get("species", ""))
        count = int(signal.value.get("count", 0))

        if species == self.x_species:
            self._current_x = count
        elif species == self.y_species:
            self._current_y = count

    def advance_to(self, t: float) -> None:
        # Record current point
        self._x_values.append(self._current_x)
        self._y_values.append(self._current_y)

        # Trim if over limit
        if len(self._x_values) > self.max_points:
            self._x_values = self._x_values[-self.max_points:]
            self._y_values = self._y_values[-self.max_points:]

        source = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "phase_point": BioSignal(
                source=source,
                name="phase_point",
                value={
                    "t": float(t),
                    "x_species": self.x_species,
                    "y_species": self.y_species,
                    "x": int(self._current_x),
                    "y": int(self._current_y),
                },
                time=float(t),
                metadata=SignalMetadata(description="Phase space current point", kind="state"),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional["VisualSpec"]:
        """Generate SVG phase space plot."""
        if len(self._x_values) < 2:
            return None

        svg = self._generate_phase_svg()
        svg_b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        data_uri = f"data:image/svg+xml;base64,{svg_b64}"

        return {
            "render": "image",
            "data": {
                "src": data_uri,
                "alt": f"Phase space: {self.x_species} vs {self.y_species}",
                "width": 500,
                "height": 400,
            },
        }

    def _generate_phase_svg(self) -> str:
        """Generate SVG phase space plot."""
        w, h = 500, 400
        margin = {"top": 30, "right": 30, "bottom": 50, "left": 60}
        plot_w = w - margin["left"] - margin["right"]
        plot_h = h - margin["top"] - margin["bottom"]

        # Determine ranges
        x_min = min(self._x_values) if self._x_values else 0
        x_max = max(self._x_values) if self._x_values else 1
        y_min = min(self._y_values) if self._y_values else 0
        y_max = max(self._y_values) if self._y_values else 1

        # Add padding
        x_range = x_max - x_min if x_max > x_min else 1
        y_range = y_max - y_min if y_max > y_min else 1
        x_min -= x_range * 0.05
        x_max += x_range * 0.05
        y_min -= y_range * 0.05
        y_max += y_range * 0.05

        def x_scale(v: float) -> float:
            return margin["left"] + ((v - x_min) / (x_max - x_min)) * plot_w

        def y_scale(v: float) -> float:
            return margin["top"] + (1 - (v - y_min) / (y_max - y_min)) * plot_h

        lines = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
            '  <style>',
            '    .axis { stroke: #333; stroke-width: 1; }',
            '    .label { font-family: sans-serif; font-size: 11px; fill: #333; }',
            '    .title { font-family: sans-serif; font-size: 13px; fill: #333; font-weight: bold; }',
            '    .trajectory { fill: none; stroke: #2563eb; stroke-width: 1.5; opacity: 0.7; }',
            '    .start { fill: #22c55e; }',
            '    .end { fill: #ef4444; }',
            '  </style>',
            f'  <rect width="{w}" height="{h}" fill="white"/>',
        ]

        # Draw axes
        x0 = margin["left"]
        x1 = margin["left"] + plot_w
        y0 = margin["top"]
        y1 = margin["top"] + plot_h

        lines.append(f'  <line class="axis" x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}"/>')
        lines.append(f'  <line class="axis" x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}"/>')

        # X-axis label
        lines.append(f'  <text class="label" x="{x0 + plot_w/2}" y="{h - 10}" text-anchor="middle">{self.x_species}</text>')

        # Y-axis label
        lines.append(f'  <text class="label" x="15" y="{y0 + plot_h/2}" text-anchor="middle" transform="rotate(-90, 15, {y0 + plot_h/2})">{self.y_species}</text>')

        # Title
        lines.append(f'  <text class="title" x="{w/2}" y="18" text-anchor="middle">Phase Space</text>')

        # Draw trajectory
        if len(self._x_values) >= 2:
            points = []
            for x, y in zip(self._x_values, self._y_values):
                px = x_scale(x)
                py = y_scale(y)
                points.append(f"{px},{py}")

            path_d = "M " + " L ".join(points)
            lines.append(f'  <path class="trajectory" d="{path_d}"/>')

            # Start point (green)
            sx, sy = x_scale(self._x_values[0]), y_scale(self._y_values[0])
            lines.append(f'  <circle class="start" cx="{sx}" cy="{sy}" r="5"/>')

            # End point (red)
            ex, ey = x_scale(self._x_values[-1]), y_scale(self._y_values[-1])
            lines.append(f'  <circle class="end" cx="{ex}" cy="{ey}" r="5"/>')

        lines.append('</svg>')
        return "\n".join(lines)

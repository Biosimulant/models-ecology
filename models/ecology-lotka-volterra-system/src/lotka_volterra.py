# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Canonical deterministic Lotka-Volterra predator-prey system."""
from __future__ import annotations

import base64
import math
from typing import Any, Dict, List, Optional, Set, TYPE_CHECKING

from biosim import BioModule
from biosim.signals import BioSignal, SignalMetadata

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


class LotkaVolterraSystem(BioModule):
    """Deterministic Lotka-Volterra predator-prey system with RK4 integration."""

    def __init__(
        self,
        alpha: float = 1.1,
        beta: float = 0.4,
        gamma: float = 0.4,
        delta: float = 0.1,
        prey_initial: float = 10.0,
        predator_initial: float = 5.0,
        prey_name: str = "Rabbits",
        predator_name: str = "Foxes",
        min_dt: float = 0.1,
    ) -> None:
        if min_dt <= 0:
            raise ValueError("min_dt must be positive")
        for name, value in {
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "delta": delta,
            "prey_initial": prey_initial,
            "predator_initial": predator_initial,
        }.items():
            if value < 0:
                raise ValueError(f"{name} must be non-negative")

        self.min_dt = float(min_dt)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.delta = float(delta)
        self.prey_initial = float(prey_initial)
        self.predator_initial = float(predator_initial)
        self.prey_name = prey_name
        self.predator_name = predator_name

        self._epsilon = 1e-9
        self._max_history_points = 10000

        self._time = 0.0
        self._prey = self.prey_initial
        self._predator = self.predator_initial
        self._history: List[Dict[str, float]] = []
        self._initial_invariant = self._invariant(self._prey, self._predator)
        self._prey_extinction_time: Optional[float] = 0.0 if self._prey <= self._epsilon else None
        self._predator_extinction_time: Optional[float] = 0.0 if self._predator <= self._epsilon else None
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return set()

    def outputs(self) -> Set[str]:
        return {"prey_population_state", "predator_population_state"}

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        self.reset()

    def reset(self) -> None:
        self._time = 0.0
        self._prey = self.prey_initial
        self._predator = self.predator_initial
        self._history = []
        self._initial_invariant = self._invariant(self._prey, self._predator)
        self._prey_extinction_time = 0.0 if self._prey <= self._epsilon else None
        self._predator_extinction_time = 0.0 if self._predator <= self._epsilon else None
        self._outputs = {}

    def advance_to(self, t: float) -> None:
        if t <= self._time:
            return

        current = self._time
        while current < t - 1e-12:
            h = min(self.min_dt, t - current)
            self._prey, self._predator = self._rk4_step(self._prey, self._predator, h)
            current += h
            self._record_state(current)

        self._time = current
        self._publish_outputs(self._time)

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def get_state(self) -> Dict[str, Any]:
        return {
            "time": self._time,
            "prey": self._prey,
            "predator": self._predator,
        }

    def visualize(self) -> Optional["VisualSpec" | List["VisualSpec"]]:
        if not self._history:
            return None

        return [
            self._population_timeseries_visual(),
            self._phase_portrait_visual(),
            self._summary_table_visual(),
            self._audit_visual(),
        ]

    def _rhs(self, prey: float, predator: float) -> tuple[float, float]:
        dprey = self.alpha * prey - self.beta * prey * predator
        dpredator = self.delta * prey * predator - self.gamma * predator
        return dprey, dpredator

    def _rk4_step(self, prey: float, predator: float, h: float) -> tuple[float, float]:
        k1x, k1y = self._rhs(prey, predator)
        k2x, k2y = self._rhs(prey + 0.5 * h * k1x, predator + 0.5 * h * k1y)
        k3x, k3y = self._rhs(prey + 0.5 * h * k2x, predator + 0.5 * h * k2y)
        k4x, k4y = self._rhs(prey + h * k3x, predator + h * k3y)

        next_prey = prey + (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        next_predator = predator + (h / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        return max(0.0, next_prey), max(0.0, next_predator)

    def _record_state(self, t: float) -> None:
        invariant = self._invariant(self._prey, self._predator)
        drift = invariant - self._initial_invariant
        self._history.append(
            {
                "t": float(t),
                "prey": float(self._prey),
                "predator": float(self._predator),
                "invariant": float(invariant),
                "drift": float(drift),
            }
        )
        if len(self._history) > self._max_history_points:
            self._history = self._history[-self._max_history_points :]

        if self._prey_extinction_time is None and self._prey <= self._epsilon:
            self._prey_extinction_time = float(t)
        if self._predator_extinction_time is None and self._predator <= self._epsilon:
            self._predator_extinction_time = float(t)

    def _publish_outputs(self, t: float) -> None:
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "prey_population_state": BioSignal(
                source=source_name,
                name="prey_population_state",
                value={"species": self.prey_name, "count": float(self._prey), "t": float(t)},
                time=float(t),
                metadata=SignalMetadata(
                    units=None,
                    description="Prey population state",
                    kind="state",
                ),
            ),
            "predator_population_state": BioSignal(
                source=source_name,
                name="predator_population_state",
                value={"species": self.predator_name, "count": float(self._predator), "t": float(t)},
                time=float(t),
                metadata=SignalMetadata(
                    units=None,
                    description="Predator population state",
                    kind="state",
                ),
            ),
        }

    def _population_timeseries_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Continuous prey and predator populations generated by the canonical Lotka-Volterra system.",
            "data": {
                "title": "Population Trajectories",
                "series": [
                    {
                        "name": self.prey_name,
                        "points": [[point["t"], point["prey"]] for point in self._history],
                    },
                    {
                        "name": self.predator_name,
                        "points": [[point["t"], point["predator"]] for point in self._history],
                    },
                ],
            },
        }

    def _phase_portrait_visual(self) -> "VisualSpec":
        svg = self._build_phase_svg()
        svg_b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
        return {
            "render": "image",
            "description": "Phase portrait with trajectory, prey and predator nullclines, and the non-trivial equilibrium.",
            "data": {
                "title": "Phase Portrait",
                "src": f"data:image/svg+xml;base64,{svg_b64}",
                "alt": "Lotka-Volterra phase portrait",
                "width": 560,
                "height": 420,
            },
        }

    def _summary_table_visual(self) -> "VisualSpec":
        prey_values = [point["prey"] for point in self._history]
        predator_values = [point["predator"] for point in self._history]
        equilibrium_prey = self.gamma / self.delta if self.delta > 0 else float("inf")
        equilibrium_predator = self.alpha / self.beta if self.beta > 0 else float("inf")
        period = self._estimate_period()
        return {
            "render": "table",
            "description": "Lotka-Volterra parameter summary, equilibrium, extrema, and extinction diagnostics.",
            "data": {
                "title": "Lotka-Volterra Summary",
                "columns": ["Metric", "Value"],
                "rows": [
                    ["alpha", f"{self.alpha:.6g}"],
                    ["beta", f"{self.beta:.6g}"],
                    ["gamma", f"{self.gamma:.6g}"],
                    ["delta", f"{self.delta:.6g}"],
                    ["Prey initial", f"{self.prey_initial:.6g}"],
                    ["Predator initial", f"{self.predator_initial:.6g}"],
                    ["Prey equilibrium X*", f"{equilibrium_prey:.6g}"],
                    ["Predator equilibrium Y*", f"{equilibrium_predator:.6g}"],
                    ["Prey final", f"{self._prey:.6g}"],
                    ["Predator final", f"{self._predator:.6g}"],
                    ["Prey min / max", f"{min(prey_values):.6g} / {max(prey_values):.6g}"],
                    ["Predator min / max", f"{min(predator_values):.6g} / {max(predator_values):.6g}"],
                    ["Estimated period", "n/a" if period is None else f"{period:.6g}"],
                    ["Prey extinction", self._format_extinction(self._prey_extinction_time)],
                    ["Predator extinction", self._format_extinction(self._predator_extinction_time)],
                ],
            },
        }

    def _audit_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Audit plot for the Lotka-Volterra conserved quantity. The drift should remain small for a faithful RK4 integration.",
            "data": {
                "title": "Invariant Drift Audit",
                "series": [
                    {
                        "name": "Invariant",
                        "points": [[point["t"], point["invariant"]] for point in self._history],
                    },
                    {
                        "name": "Drift from initial",
                        "points": [[point["t"], point["drift"]] for point in self._history],
                    },
                ],
            },
        }

    def _build_phase_svg(self) -> str:
        width, height = 560, 420
        margin_left, margin_right, margin_top, margin_bottom = 70, 24, 28, 48
        plot_width = width - margin_left - margin_right
        plot_height = height - margin_top - margin_bottom
        prey_values = [point["prey"] for point in self._history]
        predator_values = [point["predator"] for point in self._history]
        equilibrium_prey = self.gamma / self.delta if self.delta > 0 else 0.0
        equilibrium_predator = self.alpha / self.beta if self.beta > 0 else 0.0
        x_max = max(prey_values + [equilibrium_prey]) * 1.1 or 1.0
        y_max = max(predator_values + [equilibrium_predator]) * 1.1 or 1.0
        x_min = 0.0
        y_min = 0.0

        def scale_x(value: float) -> float:
            return margin_left + ((value - x_min) / max(x_max - x_min, self._epsilon)) * plot_width

        def scale_y(value: float) -> float:
            return margin_top + (1.0 - ((value - y_min) / max(y_max - y_min, self._epsilon))) * plot_height

        trajectory = " ".join(
            (
                ("M" if index == 0 else "L")
                + f" {scale_x(point['prey']):.3f},{scale_y(point['predator']):.3f}"
            )
            for index, point in enumerate(self._history)
        )
        x_eq = scale_x(equilibrium_prey)
        y_eq = scale_y(equilibrium_predator)
        y_nullcline = scale_y(equilibrium_predator)
        x_nullcline = scale_x(equilibrium_prey)
        prey_start = self._history[0]["prey"]
        predator_start = self._history[0]["predator"]
        prey_end = self._history[-1]["prey"]
        predator_end = self._history[-1]["predator"]

        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                "<style>",
                ".axis { stroke: #334155; stroke-width: 1.5; }",
                ".tick { stroke: #cbd5e1; stroke-width: 1; }",
                ".label { fill: #cbd5e1; font-size: 12px; font-family: sans-serif; }",
                ".title { fill: #f8fafc; font-size: 14px; font-family: sans-serif; font-weight: bold; }",
                ".trajectory { fill: none; stroke: #38bdf8; stroke-width: 2; }",
                ".prey-nullcline { stroke: #f97316; stroke-width: 1.5; stroke-dasharray: 6 4; }",
                ".predator-nullcline { stroke: #22c55e; stroke-width: 1.5; stroke-dasharray: 6 4; }",
                ".equilibrium { fill: #facc15; stroke: #111827; stroke-width: 1; }",
                ".start { fill: #22c55e; }",
                ".end { fill: #ef4444; }",
                "</style>",
                f'<rect width="{width}" height="{height}" fill="#0f172a" rx="12"/>',
                f'<text class="title" x="{width / 2:.1f}" y="18" text-anchor="middle">Lotka-Volterra Phase Portrait</text>',
                f'<line class="axis" x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" y2="{margin_top + plot_height}"/>',
                f'<line class="axis" x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}"/>',
                f'<line class="prey-nullcline" x1="{margin_left}" y1="{y_nullcline:.3f}" x2="{margin_left + plot_width}" y2="{y_nullcline:.3f}"/>',
                f'<line class="predator-nullcline" x1="{x_nullcline:.3f}" y1="{margin_top}" x2="{x_nullcline:.3f}" y2="{margin_top + plot_height}"/>',
                f'<path class="trajectory" d="{trajectory}"/>',
                f'<circle class="equilibrium" cx="{x_eq:.3f}" cy="{y_eq:.3f}" r="5"/>',
                f'<circle class="start" cx="{scale_x(prey_start):.3f}" cy="{scale_y(predator_start):.3f}" r="4"/>',
                f'<circle class="end" cx="{scale_x(prey_end):.3f}" cy="{scale_y(predator_end):.3f}" r="4"/>',
                f'<text class="label" x="{margin_left + plot_width / 2:.1f}" y="{height - 14}" text-anchor="middle">{self.prey_name}</text>',
                f'<text class="label" x="18" y="{margin_top + plot_height / 2:.1f}" text-anchor="middle" transform="rotate(-90, 18, {margin_top + plot_height / 2:.1f})">{self.predator_name}</text>',
                f'<text class="label" x="{margin_left + 8}" y="{y_nullcline - 8:.3f}">dX/dt = 0</text>',
                f'<text class="label" x="{x_nullcline + 8:.3f}" y="{margin_top + 16}">dY/dt = 0</text>',
                f'<text class="label" x="{x_eq + 8:.3f}" y="{y_eq - 8:.3f}">equilibrium</text>',
                "</svg>",
            ]
        )

    def _invariant(self, prey: float, predator: float) -> float:
        prey_clip = max(prey, self._epsilon)
        predator_clip = max(predator, self._epsilon)
        return (
            self.delta * prey_clip
            - self.gamma * math.log(prey_clip)
            + self.beta * predator_clip
            - self.alpha * math.log(predator_clip)
        )

    def _estimate_period(self) -> Optional[float]:
        maxima: List[float] = []
        for prev_point, point, next_point in zip(self._history, self._history[1:], self._history[2:]):
            if point["prey"] > prev_point["prey"] and point["prey"] >= next_point["prey"]:
                maxima.append(point["t"])
        if len(maxima) < 2:
            return None
        intervals = [b - a for a, b in zip(maxima, maxima[1:])]
        return sum(intervals) / len(intervals)

    def _format_extinction(self, extinction_time: Optional[float]) -> str:
        if extinction_time is None:
            return "no"
        return f"yes at t={extinction_time:.6g}"

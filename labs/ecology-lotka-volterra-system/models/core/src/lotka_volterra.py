# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Canonical deterministic Lotka-Volterra predator-prey system."""
from __future__ import annotations

import base64
import math
from typing import Any, Dict, List, Optional, TYPE_CHECKING

from biosim import StatefulBioModule
from biosim.signals import BioSignal, SignalSpec, coerce_float, scalar_or_record_input

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim.visuals import VisualSpec


class LotkaVolterraSystem(StatefulBioModule):
    """Deterministic Lotka-Volterra predator-prey system with RK4 integration."""

    def __init__(
        self,
        alpha: float = 1.1,
        beta: float = 0.4,
        gamma: float = 0.4,
        delta: float = 0.1,
        prey_initial: float = 10.0,
        predator_initial: float = 5.0,
        prey_name: str = "Prey",
        predator_name: str = "Predator",
        integration_step: float = 0.1,
    ) -> None:
        if integration_step <= 0:
            raise ValueError("integration_step must be positive")
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

        super().__init__(
            integration_step=integration_step,
            max_history_points=10000,
            publish_on_zero_window=False,
        )
        self.integration_step = float(integration_step)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.delta = float(delta)
        self.prey_initial = float(prey_initial)
        self.predator_initial = float(predator_initial)
        self.prey_name = prey_name or "Prey"
        self.predator_name = predator_name or "Predator"

        self._epsilon = 1e-9
        self._prey = self.prey_initial
        self._predator = self.predator_initial
        self._initial_invariant = self._invariant(self._prey, self._predator)
        self._prey_extinction_time: Optional[float] = 0.0 if self._prey <= self._epsilon else None
        self._predator_extinction_time: Optional[float] = 0.0 if self._predator <= self._epsilon else None

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "prey_initial_population": scalar_or_record_input("count", "Initial prey population count."),
            "predator_initial_population": scalar_or_record_input("count", "Initial predator population count."),
            "prey_growth_rate": scalar_or_record_input("1/day", "Intrinsic prey growth rate alpha."),
            "predation_rate": scalar_or_record_input("1/(count*day)", "Mass-action predation coefficient beta."),
            "predator_mortality_rate": scalar_or_record_input("1/day", "Predator mortality rate gamma."),
            "predator_reproduction_rate": scalar_or_record_input(
                "1/(count*day)", "Predator reproduction coefficient delta."
            ),
        }

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            "prey_population_state": SignalSpec.record(
                schema={"role": "str", "label": "str", "count": "float", "t": "float"},
                emitted_unit="count",
                description="Current prey population count.",
            ),
            "predator_population_state": SignalSpec.record(
                schema={"role": "str", "label": "str", "count": "float", "t": "float"},
                emitted_unit="count",
                description="Current predator population count.",
            ),
            "visualisation_payload": SignalSpec.record(
                schema={"payload": "json"},
                description="Internal history payload for the sibling visualisation model.",
            ),
        }

    def reset_state(self) -> None:
        self._prey = self.prey_initial
        self._predator = self.predator_initial
        self._initial_invariant = self._invariant(self._prey, self._predator)
        self._prey_extinction_time = 0.0 if self._prey <= self._epsilon else None
        self._predator_extinction_time = 0.0 if self._predator <= self._epsilon else None

    def _input_number(self, name: str) -> float | None:
        signal = self._input_overrides.get(name)
        if signal is None:
            return None
        return coerce_float(signal)

    def apply_overrides(self, *, reset_initial_state: bool) -> None:
        for input_name, attr_name in {
            "prey_growth_rate": "alpha",
            "predation_rate": "beta",
            "predator_mortality_rate": "gamma",
            "predator_reproduction_rate": "delta",
        }.items():
            value = self._input_number(input_name)
            if value is not None and value >= 0.0:
                setattr(self, attr_name, value)

        prey_initial = self._input_number("prey_initial_population")
        predator_initial = self._input_number("predator_initial_population")
        if prey_initial is not None and prey_initial >= 0.0:
            self.prey_initial = prey_initial
            if reset_initial_state:
                self._prey = prey_initial
        if predator_initial is not None and predator_initial >= 0.0:
            self.predator_initial = predator_initial
            if reset_initial_state:
                self._predator = predator_initial

    def get_state(self) -> Dict[str, Any]:
        return {
            "time": self._time,
            "prey": self._prey,
            "predator": self._predator,
        }

    def visualize(self) -> Optional["VisualSpec" | List["VisualSpec"]]:
        return None

    def _rhs(self, prey: float, predator: float) -> tuple[float, float]:
        dprey = self.alpha * prey - self.beta * prey * predator
        dpredator = self.delta * prey * predator - self.gamma * predator
        return dprey, dpredator

    def _derivatives(self, prey: float, predator: float) -> tuple[float, float]:
        return self._rhs(prey, predator)

    def _rk4_step(self, prey: float, predator: float, h: float) -> tuple[float, float]:
        k1x, k1y = self._rhs(prey, predator)
        k2x, k2y = self._rhs(prey + 0.5 * h * k1x, predator + 0.5 * h * k1y)
        k3x, k3y = self._rhs(prey + 0.5 * h * k2x, predator + 0.5 * h * k2y)
        k4x, k4y = self._rhs(prey + h * k3x, predator + h * k3y)

        next_prey = prey + (h / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        next_predator = predator + (h / 6.0) * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        return max(0.0, next_prey), max(0.0, next_predator)

    def step(self, h: float) -> None:
        self._prey, self._predator = self._rk4_step(self._prey, self._predator, h)

    def record_state(self, t: float) -> None:
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

        if self._prey_extinction_time is None and self._prey <= self._epsilon:
            self._prey_extinction_time = float(t)
        if self._predator_extinction_time is None and self._predator <= self._epsilon:
            self._predator_extinction_time = float(t)

    def output_payload(self, t: float) -> dict[str, Any]:
        return {
            "prey_population_state": {
                "role": "prey",
                "label": self.prey_name,
                "count": float(self._prey),
                "t": float(t),
            },
            "predator_population_state": {
                "role": "predator",
                "label": self.predator_name,
                "count": float(self._predator),
                "t": float(t),
            },
            "visualisation_payload": {"payload": self._visualisation_payload()},
        }

    def _visualisation_payload(self) -> Dict[str, Any]:
        return {
            "parameters": {
                "alpha": self.alpha,
                "beta": self.beta,
                "gamma": self.gamma,
                "delta": self.delta,
                "prey_initial": self.prey_initial,
                "predator_initial": self.predator_initial,
                "prey_name": self.prey_name,
                "predator_name": self.predator_name,
            },
            "prey_extinction_time": self._prey_extinction_time,
            "predator_extinction_time": self._predator_extinction_time,
            "history": list(self._history),
        }

    def _population_timeseries_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Canonical Lotka-Volterra prey and predator population counts.",
            "data": {
                "title": "Population Trajectories",
                "x_unit": "day",
                "y_unit": "count",
                "series": [
                    {"name": self.prey_name, "points": [[point["t"], point["prey"]] for point in self._history]},
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
            "description": "Phase portrait with trajectory, nullclines, and the non-trivial equilibrium.",
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
            "description": "Lotka-Volterra parameters, units, equilibrium, extrema, and extinction diagnostics.",
            "data": {
                "title": "Lotka-Volterra Summary",
                "columns": ["Metric", "Value", "Unit"],
                "rows": [
                    ["prey growth rate alpha", f"{self.alpha:.6g}", "1/day"],
                    ["predation rate beta", f"{self.beta:.6g}", "1/(count*day)"],
                    ["predator mortality gamma", f"{self.gamma:.6g}", "1/day"],
                    ["predator reproduction delta", f"{self.delta:.6g}", "1/(count*day)"],
                    ["prey initial", f"{self.prey_initial:.6g}", "count"],
                    ["predator initial", f"{self.predator_initial:.6g}", "count"],
                    ["prey equilibrium N*", f"{equilibrium_prey:.6g}", "count"],
                    ["predator equilibrium P*", f"{equilibrium_predator:.6g}", "count"],
                    ["prey final", f"{self._prey:.6g}", "count"],
                    ["predator final", f"{self._predator:.6g}", "count"],
                    ["prey min / max", f"{min(prey_values):.6g} / {max(prey_values):.6g}", "count"],
                    ["predator min / max", f"{min(predator_values):.6g} / {max(predator_values):.6g}", "count"],
                    ["estimated period", "n/a" if period is None else f"{period:.6g}", "day"],
                    ["prey extinction", self._format_extinction(self._prey_extinction_time), ""],
                    ["predator extinction", self._format_extinction(self._predator_extinction_time), ""],
                ],
            },
        }

    def _audit_visual(self) -> "VisualSpec":
        return {
            "render": "timeseries",
            "description": "Conserved-quantity audit. Drift should remain small for a faithful RK4 Lotka-Volterra integration.",
            "data": {
                "title": "Invariant Drift Audit",
                "series": [
                    {"name": "Invariant", "points": [[point["t"], point["invariant"]] for point in self._history]},
                    {"name": "Drift from initial", "points": [[point["t"], point["drift"]] for point in self._history]},
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

        def scale_x(value: float) -> float:
            return margin_left + (value / max(x_max, self._epsilon)) * plot_width

        def scale_y(value: float) -> float:
            return margin_top + (1.0 - value / max(y_max, self._epsilon)) * plot_height

        trajectory = " ".join(
            (("M" if index == 0 else "L") + f" {scale_x(point['prey']):.3f},{scale_y(point['predator']):.3f}")
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
                f'<text class="label" x="{margin_left + plot_width / 2:.1f}" y="{height - 14}" text-anchor="middle">{self.prey_name} count</text>',
                f'<text class="label" x="18" y="{margin_top + plot_height / 2:.1f}" text-anchor="middle" transform="rotate(-90, 18, {margin_top + plot_height / 2:.1f})">{self.predator_name} count</text>',
                f'<text class="label" x="{margin_left + 8}" y="{y_nullcline - 8:.3f}">dN/dt = 0</text>',
                f'<text class="label" x="{x_nullcline + 8:.3f}" y="{margin_top + 16}">dP/dt = 0</text>',
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

    @staticmethod
    def _format_extinction(extinction_time: Optional[float]) -> str:
        if extinction_time is None:
            return "no"
        return f"yes at t={extinction_time:.6g}"

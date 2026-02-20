# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Auto-generated SBML BioModule wrapper for Poliquin2013 - Energy Deregulations in Parkinson's Disease.

Source: biomodels_ebi:MODEL1410060000
Original: https://www.ebi.ac.uk/biomodels/MODEL1410060000
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim import BioWorld

import biosim
from biosim.signals import BioSignal, SignalMetadata


class SbmlPoliquin2013EnergyDeregulationsInParkinsonSDisease(biosim.BioModule):
    """BioModule wrapper for SBML model: Poliquin2013 - Energy Deregulations in Parkinson's Disease."""

    def __init__(self, model_path: str = "data/MODEL1410060000.xml", min_dt: float = 0.01) -> None:
        self.min_dt = min_dt
        self._model_path = Path(__file__).parent.parent / model_path
        self._t = 0.0
        self._rr = None
        self._species_ids: list[str] = []
        self._stub: bool = False
        self._outputs: Dict[str, BioSignal] = {}

    @staticmethod
    def _is_placeholder_file(path: Path) -> bool:
        try:
            head = path.open("rb").read(160)
        except OSError:
            return True
        try:
            text = head.decode("utf-8", errors="ignore")
        except Exception:
            return False
        return "Placeholder:" in text

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        if self._is_placeholder_file(self._model_path):
            # "No downloads" mode: keep the module runnable without external simulators.
            self._stub = True
            self._rr = None
            self._species_ids = []
            self._t = 0.0
            return

        try:
            import tellurium as te
        except Exception:
            self._stub = True
            self._rr = None
            self._species_ids = []
            self._t = 0.0
            return

        self._rr = te.loadSBMLModel(str(self._model_path))
        self._species_ids = list(self._rr.getFloatingSpeciesIds())
        self._t = 0.0

    def reset(self) -> None:
        if self._rr is not None:
            self._rr.reset()
        self._t = 0.0
        self._outputs = {}

    def inputs(self) -> Set[str]:
        return set()

    def outputs(self) -> Set[str]:
        return {"state"}

    def advance_to(self, t: float) -> None:
        if self._rr is None and not self._stub:
            self.setup()
        if self._stub:
            self._t = t
            source_name = getattr(self, "_world_name", self.__class__.__name__)
            self._outputs = {
                "state": BioSignal(
                    source=source_name,
                    name="state",
                    value={"t": t, "placeholder": True},
                    time=t,
                    metadata=SignalMetadata(
                        description="Placeholder SBML state (no downloads / simulator unavailable)",
                        kind="state",
                    ),
                )
            }
            return

        if t > self._t:
            self._rr.integrate(self._t, t)
            self._t = t
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        concentrations = {}
        for sid in self._species_ids:
            try:
                concentrations[sid] = float(self._rr[sid])
            except Exception:
                concentrations[sid] = 0.0
        self._outputs = {
            "state": BioSignal(
                source=source_name,
                name="state",
                value=concentrations,
                time=t,
                metadata=SignalMetadata(
                    units="mM",
                    description="Species concentrations",
                    kind="state",
                ),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

# Canonical alias for stable entrypoint naming.
Poliquin2013EnergyDeregulationsInParkinsonSModel1410060000Model = SbmlPoliquin2013EnergyDeregulationsInParkinsonSDisease

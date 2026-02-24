# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Auto-generated SBML BioModule wrapper for Geci2022.

Source: biomodels_ebi:MODEL2301120001
Original: https://www.ebi.ac.uk/biomodels/MODEL2301120001
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Set, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from biosim import BioWorld

import biosim
from biosim.signals import BioSignal, SignalMetadata

import logging

logger = logging.getLogger(__name__)

class SbmlGeci2022(biosim.BioModule):
    """BioModule wrapper for SBML model: Geci2022."""

    def __init__(self, model_path: str = "data/MODEL2301120001.xml", min_dt: float = 0.01) -> None:
        self.min_dt = min_dt
        self._model_path = Path(__file__).parent.parent / model_path
        self._t = 0.0
        self._rr = None
        self._species_ids: list[str] = []
        self._outputs: Dict[str, BioSignal] = {}

    def setup(self, config: Optional[Dict[str, Any]] = None) -> None:
        import tellurium as te

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
        if self._rr is None:
            self.setup()

        if t > self._t:
            self._rr.simulate(self._t, t)
            self._t = t

        source_name = getattr(self, "_world_name", self.__class__.__name__)
        concentrations = {}
        for sid in self._species_ids:
            try:
                concentrations[sid] = float(self._rr[sid])
            except (KeyError, ValueError, TypeError):  # narrowed from bare Exception
                logger.warning("Failed to read species %s, defaulting to 0.0", sid)
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
Geci2022Model2301120001Model = SbmlGeci2022

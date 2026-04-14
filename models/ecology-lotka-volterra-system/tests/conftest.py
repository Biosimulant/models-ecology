from __future__ import annotations

import sys
from pathlib import Path

import pytest

_MODEL_DIR = Path(__file__).resolve().parents[1]
_MONOREPO_ROOT = _MODEL_DIR.parents[3]
_BSIM_SRC = _MONOREPO_ROOT / "bsim-active" / "biosim" / "src"


@pytest.fixture(scope="session", autouse=True)
def _paths():
    for path in (str(_MODEL_DIR), str(_BSIM_SRC)):
        if path not in sys.path:
            sys.path.insert(0, path)


@pytest.fixture(scope="session")
def biosim(_paths):
    import biosim as _bsim

    return _bsim

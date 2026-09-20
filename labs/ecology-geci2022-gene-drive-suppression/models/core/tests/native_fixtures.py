"""Retrieve exact native Julia reference bytes when binary fixtures are external.

The reference data are committed outputs of the unchanged upstream Julia model,
not newly generated expectations from the Python implementation under test.
"""
import hashlib
import json
from pathlib import Path
import tempfile
import urllib.request

FIXTURES = Path(__file__).resolve().parents[3] / "verification"
REFERENCE_COMMIT = "ad82fc61a6a67843ebca8911a4bf6ffb8627e407"
REFERENCE_URL = (
    "https://raw.githubusercontent.com/Biosimulant/models-ecology/"
    + REFERENCE_COMMIT
    + "/labs/ecology-geci2022-gene-drive-suppression/verification/"
)
MAX_BYTES = 2 * 1024 * 1024


def native_fixture(filename: str, *, directory: Path = FIXTURES) -> Path:
    expected = json.loads((directory / "fixture-sha256.json").read_text())
    if filename not in expected or Path(filename).name != filename:
        raise ValueError("requested native fixture is not in the pinned manifest")
    target = directory / filename
    if target.exists():
        data = target.read_bytes()
    else:
        with urllib.request.urlopen(REFERENCE_URL + filename, timeout=60) as response:
            data = response.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES or hashlib.sha256(data).hexdigest() != expected[filename]:
        raise ValueError(f"native reference integrity failure: {filename}")
    if not target.exists():
        with tempfile.NamedTemporaryFile(dir=directory, prefix=".native-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(data)
        try:
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
    return target

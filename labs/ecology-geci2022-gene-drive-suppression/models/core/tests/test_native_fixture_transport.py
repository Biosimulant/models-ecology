import hashlib
import io
import json

import pytest
import native_fixtures


def test_missing_fixture_fetches_exact_pinned_bytes_and_reuses_cache(tmp_path, monkeypatch):
    data = b"reference bytes"
    (tmp_path / "fixture-sha256.json").write_text(json.dumps({"example.tsv.gz": hashlib.sha256(data).hexdigest()}))
    urls = []

    def fetch(url, **kwargs):
        urls.append(url)
        return io.BytesIO(data)

    monkeypatch.setattr(native_fixtures.urllib.request, "urlopen", fetch)
    assert native_fixtures.native_fixture("example.tsv.gz", directory=tmp_path).read_bytes() == data
    native_fixtures.native_fixture("example.tsv.gz", directory=tmp_path)
    assert urls == [native_fixtures.REFERENCE_URL + "example.tsv.gz"]


@pytest.mark.parametrize("cached", [False, True])
def test_changed_reference_bytes_fail_instead_of_changing_expectations(tmp_path, monkeypatch, cached):
    (tmp_path / "fixture-sha256.json").write_text(json.dumps({"example.tsv.gz": "0" * 64}))
    if cached:
        (tmp_path / "example.tsv.gz").write_bytes(b"changed")
    monkeypatch.setattr(native_fixtures.urllib.request, "urlopen", lambda *args, **kwargs: io.BytesIO(b"changed"))
    with pytest.raises(ValueError, match="integrity failure"):
        native_fixtures.native_fixture("example.tsv.gz", directory=tmp_path)


def test_unlisted_reference_cannot_be_fetched(tmp_path):
    (tmp_path / "fixture-sha256.json").write_text("{}")
    with pytest.raises(ValueError, match="pinned manifest"):
        native_fixtures.native_fixture("../other", directory=tmp_path)

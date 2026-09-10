import hashlib

import pytest

from yt.utilities.answer_testing.testing_utilities import generate_hash


@pytest.mark.parametrize(
    ("data", "payload"),
    [(b"yt", b"yt"), (None, b"-1"), ({"a": b"b"}, b"ab")],
)
def test_generate_hash_in_fips_mode(monkeypatch, data, payload):
    expected = hashlib.md5(payload, usedforsecurity=False).hexdigest()
    original_md5 = hashlib.md5

    def fips_md5(data=b"", *, usedforsecurity=True):
        if usedforsecurity:
            raise ValueError("MD5 is disabled for FIPS")
        return original_md5(data)

    monkeypatch.setattr(hashlib, "md5", fips_md5)

    assert generate_hash(data) == expected

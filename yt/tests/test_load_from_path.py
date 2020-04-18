import os
from sys import version_info
from yt.convenience import _sanitize_load_args
import pytest

PY36 = version_info >= (3, 6)

if PY36:
    from pathlib import Path

@pytest.mark.skipif(not PY36, reason="requires python3.6 or higher")
class TestSanitizeLoadArgs():
    def test_sanitize_from_simple_path(self):
        # check for basic path
        p1 = os.path.join("not", "a", "real", "datafile.hdf5")
        p2 = Path(p1)
        assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

    def test_sanitize_from_two_paths(self):
        # check for more than one path
        p1 = [os.path.join("not", "a", "real", "datafile.hdf5"), os.path.join("not", "real", "either", "datafile.hdf5")]
        p2 = [Path(p) for p in p1]
        assert _sanitize_load_args(*p1) == _sanitize_load_args(*p2)

    def test_sanitize_from_user_path(self):
        # check for user "~" card expansion
        p1 = os.path.join("~", "not", "a", "real", "datafile.hdf5")
        p2 = Path(p1)
        assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

    def test_sanitize_from_wildcard_path(self):
        # check with wildcards
        p1 = os.path.join("not", "a", "real", "directory", "*.hdf5")
        p2 = Path(p1)
        assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

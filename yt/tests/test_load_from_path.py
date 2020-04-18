from yt.convenience import _sanitize_load_args
from pathlib import Path

def test_sanitize_from_simple_path():
    # check for basic path
    p1 = "not/a/real/datafile.hdf5"
    p2 = Path(p1)
    assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

def test_sanitize_from_two_paths():
    # check for more than one path
    p1 = ["not/a/real/datafile.hdf5", "not/real/either/datafile.hdf5"]
    p2 = [Path(p) for p in p1]
    assert _sanitize_load_args(*p1) == _sanitize_load_args(*p2)

def test_sanitize_from_user_path():
    # check for user "~" card expansion
    p1 = "~/not/a/real/datafile.hdf5"
    p2 = Path(p1)
    assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

def test_sanitize_from_wildcard_path():
    # check with wildcards
    p1 = "not/a/real/directory/*.hdf5"
    p2 = Path(p1)
    assert _sanitize_load_args(p1) == _sanitize_load_args(p2)

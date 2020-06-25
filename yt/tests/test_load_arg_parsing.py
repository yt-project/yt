import json
import os
from pathlib import Path

import pytest

from yt.convenience import _sanitize_load_args, load
from yt.funcs import ensure_list
from yt.utilities.exceptions import YTOutputNotIdentified


def test_sanitize_empty():
    sargs = _sanitize_load_args()
    assert sargs == []


@pytest.mark.parametrize("obj", [None, 123, 2.0, (None, 123), [None, 123]])
def test_sanitize_arbitrary_obj(obj):
    # check that (nested) objects which are neither pathlike nor iterable are unaltered
    args = ensure_list(obj)
    sargs = _sanitize_load_args(*args)

    # this is a quick hack to check that results are structuraly equivalent
    # if tuples are converted to list or vice versa it's not expected to be an issue
    assert json.dumps(sargs) == json.dumps(args)


def test_raise_filenotfound(tmpdir):
    with pytest.raises(FileNotFoundError):
        load(os.path.join(tmpdir, "non_existing_file"))


def test_raise_unidentifiedtype(tmpdir):
    p = tmpdir.mkdir("sub").join("invalid_data.txt")
    p.write("This should never be valid")
    with pytest.raises(YTOutputNotIdentified):
        load(p)


class TestSanitizePathLike:
    def test_sanitize_from_simple_path(self):
        # check for basic path
        p1 = os.path.join("not", "a", "real", "datafile.hdf5")

        sp1 = _sanitize_load_args(p1)

        assert type(sp1) is list
        assert type(sp1[0]) is str
        assert sp1 == _sanitize_load_args(Path(p1))

    def test_sanitize_from_two_paths(self):
        # check for more than one path
        p1 = [
            os.path.join("not", "a", "real", "datafile.hdf5"),
            os.path.join("not", "real", "either", "datafile.hdf5"),
        ]

        sp1 = _sanitize_load_args(*p1)
        assert sp1 == p1

        p2 = [Path(p) for p in p1]
        assert sp1 == _sanitize_load_args(*p2)

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

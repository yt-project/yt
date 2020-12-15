import difflib
import os
import tempfile
from pathlib import Path
from pprint import pprint

import pytest
from nose import SkipTest

# local
from yt.frontends.idefix.inifile_io import (
    IDEFIX_INI_SCHEMA,
    IdefixConf,
    read_idefix_inifile,
    write_idefix_inifile,
)

msg = "skipping idefix unit tests (environment variable $IDEFIX_DIR is not set)."
if "IDEFIX_DIR" not in os.environ:
    raise SkipTest(msg)
    # within pytest:
    # pytest.skip(msg, allow_module_level=True)

IDEFIX_DIR = Path(os.environ["IDEFIX_DIR"])
TEST_DIR = Path.joinpath(IDEFIX_DIR, "test")
IDEFIX_INI_FILES = list(TEST_DIR.glob("**/*.ini"))


@pytest.mark.parametrize(
    "invalid_data", ["", "SingleToken", "[InvalidSectionName", "InvalidSection]"]
)
def test_tokenizer(invalid_data):
    with pytest.raises(ValueError):
        IdefixConf.tokenize_line(invalid_data)


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_unit_read(inifile):
    data = IdefixConf(inifile)
    pprint(data)


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_read(inifile):
    data = read_idefix_inifile(inifile)
    pprint(data)


def test_oop_write():
    inifile = IDEFIX_INI_FILES[0]
    conf = IdefixConf(inifile)
    with tempfile.TemporaryFile(mode="wt") as tmpfile:
        conf.write(tmpfile)


def test_func_write():
    inifile = IDEFIX_INI_FILES[0]
    data = IdefixConf(inifile)
    with tempfile.TemporaryFile(mode="wt") as tmpfile:
        write_idefix_inifile(data, tmpfile)


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_idempotent_io(inifile):
    data0 = read_idefix_inifile(inifile)
    with tempfile.TemporaryDirectory() as tmpdir:
        save1 = Path(tmpdir) / "save1"
        save2 = Path(tmpdir) / "save2"
        write_idefix_inifile(data0, save1)
        data1 = read_idefix_inifile(save1)
        write_idefix_inifile(data1, save2)

        text1 = open(save1, "r").readlines()
        text2 = open(save2, "r").readlines()

        diff = "".join(difflib.context_diff(text1, text2))
        assert not diff


@pytest.mark.parametrize("func", [read_idefix_inifile, write_idefix_inifile])
def test_docstrings(func):
    assert IDEFIX_INI_SCHEMA in func.__doc__

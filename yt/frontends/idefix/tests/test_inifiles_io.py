import filecmp
import os
import tempfile
from pathlib import Path
from pprint import pprint

import pytest

# local
from yt.frontends.idefix.inifile_io import (
    IDEFIX_INI_SCHEMA,
    IdefixConf,
    read_idefix_inifile,
    write_idefix_inifile,
)

if "IDEFIX_DIR" not in os.environ:
    pytest.skip(
        "skipping idefix inifile unit tests "
        "(reason: environment variable $IDEFIX_DIR is not set).",
        allow_module_level=True,
    )

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


def test_unit_write():
    inifile = IDEFIX_INI_FILES[0]
    conf = IdefixConf(inifile)
    conf.write("test1.ini.tmp")


def test_write():
    inifile = IDEFIX_INI_FILES[0]
    data = IdefixConf(inifile)
    write_idefix_inifile(data, "test2.ini.tmp")


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_idempotent_io(inifile):
    data0 = read_idefix_inifile(inifile)
    with tempfile.TemporaryDirectory() as tmpdir:
        save1 = Path(tmpdir) / "save1"
        save2 = Path(tmpdir) / "save2"
        write_idefix_inifile(data0, save1)
        data1 = read_idefix_inifile(save1)
        write_idefix_inifile(data1, save2)
        assert filecmp.cmp(save1, save2)


@pytest.mark.parametrize("func", [read_idefix_inifile, write_idefix_inifile])
def test_docstrings(func):
    assert IDEFIX_INI_SCHEMA in func.__doc__

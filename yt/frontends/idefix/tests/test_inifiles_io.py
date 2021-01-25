import difflib
import os
import tempfile
from pathlib import Path

import pytest

from yt.frontends.idefix.inifile_io import IdefixConf

msg = "skipping idefix unit tests (environment variable $IDEFIX_DIR is not set)."
if "IDEFIX_DIR" not in os.environ:
    from nose import SkipTest

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
    IdefixConf(inifile)


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_oop_write(inifile):
    conf = IdefixConf(inifile)
    with tempfile.TemporaryFile(mode="wt") as tmpfile:
        conf.write(tmpfile)


@pytest.mark.parametrize("inifile", IDEFIX_INI_FILES)
def test_idempotent_io(inifile):
    data0 = IdefixConf(inifile)
    with tempfile.TemporaryDirectory() as tmpdir:
        save1 = Path(tmpdir) / "save1"
        save2 = Path(tmpdir) / "save2"
        data0.write(save1)
        data1 = IdefixConf(save1)
        data1.write(save2)

        text1 = open(save1).readlines()
        text2 = open(save2).readlines()

        diff = "".join(difflib.context_diff(text1, text2))
        assert not diff

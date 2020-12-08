import os
from pathlib import Path

from yt.frontends.idefix.api import IdefixDataset
from yt.loaders import load

dmpfile = Path(os.environ["IDEFIX_DIR"]).joinpath("test", "HD", "KHI", "dump.0001.dmp")
inifile = dmpfile.parent / "idefix.ini"


def test_load():
    ds = load(dmpfile, inifile=inifile)
    assert isinstance(ds, IdefixDataset)
    assert ds.dimensionality == 2


def test_region():
    ds = load(dmpfile, inifile=inifile)
    ds.r[:]


def test_fields():
    ds = load(dmpfile, inifile=inifile)
    expected = [("idefix", "Vc-RHO"), ("idefix", "Vc-VX1"), ("idefix", "Vc-VX2")]
    assert ds.field_list == expected

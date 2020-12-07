import os

from yt.frontends.idefix.api import IdefixDataset
from yt.loaders import load

test_file = os.path.join(os.environ["IDEFIX_DIR"], "test", "HD", "KHI", "dump.0001.dmp")


def test_load():
    ds = load(test_file)
    assert isinstance(ds, IdefixDataset)


def test_region():
    ds = load(test_file)
    ds.r[:]

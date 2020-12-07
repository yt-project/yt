import os

from yt.loaders import load

test_file = os.path.join(os.environ["IDEFIX_DIR"], "test", "HD", "KHI", "dump.0001.dmp")


def test_load():
    load(test_file)

from yt.testing import (
    fake_random_ds,
    requires_file,
    requires_module,
)
from yt.utilities.answer_testing.framework import data_dir_load

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_module("h5py")
@requires_file(g30)
def test_checksum():
    assert fake_random_ds(16).checksum == "notafile"
    assert data_dir_load(g30).checksum == "6169536e4b9f737ce3d3ad440df44c58"

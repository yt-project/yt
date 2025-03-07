import os

import pytest

import yt
from yt.config import ytcfg
from yt.utilities.parameter_file_storage import ParameterFileStore


@pytest.fixture
def store_parameter_files():
    spf_init = ytcfg.get("yt", "store_parameter_files")
    ytcfg.set("yt", "store_parameter_files", True)
    yield
    ytcfg.set("yt", "store_parameter_files", spf_init)


@pytest.fixture
def set_parameter_file(tmp_path):
    pfs_init_name = ytcfg.get("yt", "parameter_file_store")
    pfs_path = tmp_path / "param_file_store.csv"
    pfs_path.touch()
    ytcfg.set("yt", "parameter_file_store", str(pfs_path))
    yield
    ytcfg.set("yt", "parameter_file_store", pfs_init_name)


def test_on_disk_parameter_file_store(
    monkeypatch, store_parameter_files, set_parameter_file
):
    pfs = ParameterFileStore()
    # patching _register to True forces a full re-initialization on next
    # instantiation of a ParameterFileStore
    monkeypatch.setattr(pfs, "_register", True)
    monkeypatch.setattr(pfs, "_records", {})

    # to force a re-initialization, set _register to True then get a new
    # instance. This should read from disk.
    pfs = ParameterFileStore()
    pfs_path = ytcfg.get("yt", "parameter_file_store")
    assert os.path.isfile(pfs_path)  # on disk store should now exist
    db_records = pfs.read_db()  # force a read from disk
    assert len(db_records) == 0  # and it should be empty

    # check that ds load is registered on disk
    ds = yt.load_sample("IsolatedGalaxy")
    db_records = pfs.read_db()

    assert len(db_records) == 1
    hash = ds._hash()
    assert hash in db_records
    ds_rec = db_records[hash]
    assert ds_rec["bn"] == "galaxy0030"
    assert ds_rec["class_name"] == "EnzoDataset"
    assert ds_rec["module_name"] == "yt.frontends.enzo.data_structures"

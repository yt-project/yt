import os
import tempfile
from pathlib import Path
from unittest import mock

import pytest

import yt.utilities.parameter_file_storage as pfs


class _FakeDataset:
    def __init__(
        self,
        hash_value,
        basename="output_0001",
        directory="/tmp/data",
        current_time=1.25,
        unique_identifier="ctid-1",
        instantiated=10.0,
    ):
        self._hash_value = hash_value
        self.basename = basename
        self.directory = directory
        self.current_time = current_time
        self.unique_identifier = unique_identifier
        self._instantiated = instantiated

    def _hash(self):
        return self._hash_value


class FakeLoadedDataset:
    def __init__(self, path):
        self.path = path
        self._instantiated = 99.0

    def _hash(self):
        return "hash-1"


def _bare_store():
    store = object.__new__(pfs.ParameterFileStore)
    store._distributed = False
    store._processing = False
    store._owner = 0
    store._read_only = False
    store._records = {}
    return store


def test_unknown_dataset_type_and_store_initialization_modes():
    err = pfs.UnknownDatasetType("MyDataset")
    assert str(err) == "MyDataset"
    assert repr(err) == "MyDataset"

    with (
        mock.patch.object(pfs.ParameterFileStore, "_shared_state", {}),
        mock.patch.object(pfs.ParameterFileStore, "_register", True),
        mock.patch.object(pfs.ytcfg, "get", return_value=False),
    ):
        store1 = pfs.ParameterFileStore()
        store2 = pfs.ParameterFileStore()

    assert store1.__dict__ is store2.__dict__
    assert store1._read_only is True
    assert store1._records == {}

    with (
        mock.patch.object(pfs.ParameterFileStore, "_shared_state", {}),
        mock.patch.object(pfs.ParameterFileStore, "_register", True),
        mock.patch.object(pfs.ytcfg, "get", return_value=True),
        mock.patch.object(pfs.ParameterFileStore, "init_db") as init_db,
        mock.patch.object(
            pfs.ParameterFileStore, "read_db", return_value={"hash-1": {"bn": "x"}}
        ) as read_db,
    ):
        writable_store = pfs.ParameterFileStore()

    init_db.assert_called_once_with()
    read_db.assert_called_once_with()
    assert writable_store._read_only is False
    assert writable_store._records == {"hash-1": {"bn": "x"}}


def test_init_db_and_get_db_name_cover_creation_and_error_paths():
    store = _bare_store()

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "nested" / "datasets.csv"
        with mock.patch.object(store, "_get_db_name", return_value=str(db_path)):
            store.init_db()

        assert db_path.parent.is_dir()
        assert db_path.is_file()

        existing_store = _bare_store()
        existing_path = Path(tmpdir) / "existing" / "datasets.csv"
        existing_path.parent.mkdir()
        with (
            mock.patch.object(existing_store, "_get_db_name", return_value=str(existing_path)),
            mock.patch.object(pfs.os, "mkdir") as mkdir_mock,
        ):
            existing_store.init_db()
        mkdir_mock.assert_not_called()
        assert existing_path.is_file()

    with (
        mock.patch.object(store, "_get_db_name", return_value="/tmp/test/store.csv"),
        mock.patch.object(pfs.os.path, "isdir", return_value=False),
        mock.patch.object(pfs.os, "mkdir", side_effect=OSError("denied")),
        pytest.raises(pfs.NoParameterShelf),
    ):
        store.init_db()

    with mock.patch.object(pfs.ytcfg, "get", return_value="store.csv"):
        with mock.patch.object(pfs.os, "access", return_value=False):
            assert store._get_db_name() == os.path.abspath("store.csv")
        with mock.patch.object(pfs.os, "access", return_value=True):
            assert store._get_db_name().endswith(os.path.join(".yt", "store.csv"))


def test_dataset_lookup_adaptation_and_conversion_paths():
    store = _bare_store()
    record = {
        "bn": "output_0001",
        "fp": "/tmp/data",
        "tt": 1.25,
        "ctid": "ctid-1",
        "class_name": "FakeLoadedDataset",
        "last_seen": 10.0,
    }
    store._records = {"hash-1": record.copy(), "hash-2": {**record, "ctid": "ctid-2"}}

    with mock.patch.object(store, "_convert_ds", side_effect=["from-hash", "from-ctid"]) as convert:
        assert store.get_ds_hash("hash-1") == "from-hash"
        assert store.get_ds_ctid("ctid-2") == "from-ctid"
        assert store.get_ds_ctid("missing") is None

    assert convert.call_count == 2

    ds = _FakeDataset("hash-1")
    assert store._adapt_ds(ds) == {
        "bn": "output_0001",
        "fp": "/tmp/data",
        "tt": 1.25,
        "ctid": "ctid-1",
        "class_name": "_FakeDataset",
        "last_seen": 10.0,
    }

    with pytest.raises(pfs.UnknownDatasetType):
        store._convert_ds({**record, "class_name": "MissingType"})

    with (
        mock.patch.dict(pfs.output_type_registry, {"FakeLoadedDataset": FakeLoadedDataset}, clear=True),
        mock.patch.object(pfs.mylog, "info") as info_log,
        mock.patch.object(pfs.os.path, "exists", return_value=True),
    ):
        loaded = store._convert_ds(record)

    info_log.assert_called_once_with("Checking %s", "/tmp/data/output_0001")
    assert isinstance(loaded, FakeLoadedDataset)
    assert loaded.path == "/tmp/data/output_0001"
    assert store._records["hash-1"]["last_seen"] == 99.0

    with (
        mock.patch.dict(pfs.output_type_registry, {"FakeLoadedDataset": FakeLoadedDataset}, clear=True),
        mock.patch.object(pfs.os.path, "exists", return_value=False),
        pytest.raises(OSError),
    ):
        store._convert_ds(record)


def test_check_insert_wipe_flush_and_recent_behaviors():
    new_store = _bare_store()
    ds = _FakeDataset("hash-1", instantiated=12.0)
    with mock.patch.object(new_store, "insert_ds") as insert_ds:
        new_store.check_ds(ds)
    insert_ds.assert_called_once_with(ds)

    existing_store = _bare_store()
    existing_store._records = {
        "hash-1": {
            "bn": "output_0001",
            "fp": "/tmp/data",
            "tt": 1.25,
            "ctid": "ctid-1",
            "class_name": "_FakeDataset",
            "last_seen": 1.0,
        }
    }
    existing_store.check_ds(_FakeDataset("hash-1", instantiated=15.0))
    assert existing_store._records["hash-1"]["last_seen"] == 15.0

    moved_store = _bare_store()
    moved_store._records = {
        "hash-1": {
            "bn": "old",
            "fp": "/tmp/old",
            "tt": 1.25,
            "ctid": "ctid-1",
            "class_name": "_FakeDataset",
            "last_seen": 1.0,
        }
    }
    moved_ds = _FakeDataset("hash-1", basename="new", directory="/tmp/new")
    with (
        mock.patch.object(moved_store, "wipe_hash") as wipe_hash,
        mock.patch.object(moved_store, "insert_ds") as insert_ds,
    ):
        moved_store.check_ds(moved_ds)
    wipe_hash.assert_called_once_with("hash-1")
    insert_ds.assert_called_once_with(moved_ds)

    insert_store = _bare_store()
    with mock.patch.object(insert_store, "flush_db") as flush_db:
        insert_store.insert_ds(ds)
    flush_db.assert_called_once_with()
    assert insert_store._records["hash-1"]["bn"] == "output_0001"

    wipe_store = _bare_store()
    with mock.patch.object(wipe_store, "flush_db") as flush_db:
        wipe_store.wipe_hash("missing")
    flush_db.assert_not_called()

    wipe_store._records = {"hash-1": {"last_seen": 1.0}}
    with mock.patch.object(wipe_store, "flush_db") as flush_db:
        wipe_store.wipe_hash("hash-1")
    flush_db.assert_called_once_with()
    assert wipe_store._records == {}

    read_only_store = _bare_store()
    read_only_store._read_only = True
    with (
        mock.patch.object(read_only_store, "_write_out") as write_out,
        mock.patch.object(read_only_store, "read_db") as read_db,
    ):
        read_only_store.flush_db()
    write_out.assert_not_called()
    read_db.assert_not_called()

    writable_store = _bare_store()
    with (
        mock.patch.object(writable_store, "_write_out") as write_out,
        mock.patch.object(writable_store, "read_db", return_value={}) as read_db,
    ):
        writable_store.flush_db()
    write_out.assert_called_once_with()
    read_db.assert_called_once_with()

    writable_store._records = {
        "h1": {"last_seen": 1.0},
        "h2": {"last_seen": 4.0},
        "h3": {"last_seen": 2.0},
    }
    assert writable_store.get_recent(2) == [{"last_seen": 4.0}, {"last_seen": 2.0}]


def test_write_out_and_read_db_cover_csv_round_trip_and_limits():
    store = _bare_store()

    store._read_only = True
    assert store._write_out() is None
    store._read_only = False

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "datasets.csv"
        store._records = {
            "h1": {
                "bn": "a",
                "fp": "/tmp/a",
                "tt": "1.0",
                "ctid": "1",
                "class_name": "A",
                "last_seen": 2.0,
            },
            "h2": {
                "bn": "b",
                "fp": "/tmp/b",
                "tt": "2.0",
                "ctid": "2",
                "class_name": "B",
                "last_seen": 5.0,
            },
        }

        with (
            mock.patch.object(store, "_get_db_name", return_value=str(db_path)),
            mock.patch.object(pfs.ytcfg, "get", return_value=1),
        ):
            store._write_out()

        assert db_path.read_text(encoding="utf-8") == "h2,b,/tmp/b,2.0,2,B,5.0\n"
        assert store._records["h2"]["hash"] == "h2"

        db_path.write_text(
            "h1,a,/tmp/a,1.0,1,A\n"
            "h2,b,/tmp/b,2.0,2,B,5.5\n",
            encoding="utf-8",
        )

        with mock.patch.object(store, "_get_db_name", return_value=str(db_path)):
            records = store.read_db()

    assert records == {
        "h1": {
            "bn": "a",
            "fp": "/tmp/a",
            "tt": "1.0",
            "ctid": "1",
            "class_name": "A",
            "last_seen": 0.0,
        },
        "h2": {
            "bn": "b",
            "fp": "/tmp/b",
            "tt": "2.0",
            "ctid": "2",
            "class_name": "B",
            "last_seen": 5.5,
        },
    }

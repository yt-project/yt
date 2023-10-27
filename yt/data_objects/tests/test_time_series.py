import tempfile
from pathlib import Path

import pytest
from numpy.testing import assert_raises

from yt.data_objects.static_output import Dataset
from yt.data_objects.time_series import DatasetSeries
from yt.utilities.exceptions import YTUnidentifiedDataType
from yt.utilities.object_registries import output_type_registry


def test_pattern_expansion():
    file_list = [f"fake_data_file_{str(i).zfill(4)}" for i in range(10)]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        for file in file_list:
            (tmp_path / file).touch()

        pattern = tmp_path / "fake_data_file_*"
        expected = [str(tmp_path / file) for file in file_list]
        found = DatasetSeries._get_filenames_from_glob_pattern(pattern)
        assert found == expected

        found2 = DatasetSeries._get_filenames_from_glob_pattern(Path(pattern))
        assert found2 == expected


def test_no_match_pattern():
    with tempfile.TemporaryDirectory() as tmpdir:
        pattern = Path(tmpdir).joinpath("fake_data_file_*")
        assert_raises(
            FileNotFoundError, DatasetSeries._get_filenames_from_glob_pattern, pattern
        )


@pytest.fixture
def FakeDataset():
    i = 0

    class __FakeDataset(Dataset):
        """A minimal loadable fake dataset subclass"""

        def __init__(self, *args, **kwargs):
            nonlocal i
            super().__init__(*args, **kwargs)
            self.current_time = i
            self.current_opposite_time = -i
            i += 1

        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

        def _parse_parameter_file(self):
            return

        def _set_code_unit_attributes(self):
            return

        def set_code_units(self):
            self.current_time = 0
            return

        def _hash(self):
            return

        def _setup_classes(self):
            return

    try:
        yield __FakeDataset
    finally:
        output_type_registry.pop("__FakeDataset")


@pytest.fixture
def fake_datasets():
    file_list = [f"fake_data_file_{str(i).zfill(4)}" for i in range(10)]
    with tempfile.TemporaryDirectory() as tmpdir:
        pfile_list = [Path(tmpdir) / file for file in file_list]
        sfile_list = [str(file) for file in pfile_list]
        for file in pfile_list:
            file.touch()
        pattern = Path(tmpdir) / "fake_data_file_*"

        yield file_list, pfile_list, sfile_list, pattern


def test_init_fake_dataseries(fake_datasets):
    file_list, pfile_list, sfile_list, pattern = fake_datasets

    # init from str pattern
    ts = DatasetSeries(pattern)
    assert ts._pre_outputs == sfile_list

    # init from Path pattern
    ppattern = Path(pattern)
    ts = DatasetSeries(ppattern)
    assert ts._pre_outputs == sfile_list

    # init form str list
    ts = DatasetSeries(sfile_list)
    assert ts._pre_outputs == sfile_list

    # init form Path list
    ts = DatasetSeries(pfile_list)
    assert ts._pre_outputs == pfile_list

    # rejected input type (str repr of a list) "[file1, file2, ...]"
    assert_raises(FileNotFoundError, DatasetSeries, str(file_list))

    # finally, check that ts[0] fails to actually load
    assert_raises(YTUnidentifiedDataType, ts.__getitem__, 0)


def test_init_fake_dataseries2(FakeDataset, fake_datasets):
    _file_list, _pfile_list, _sfile_list, pattern = fake_datasets
    ds = DatasetSeries(pattern)[0]
    assert isinstance(ds, FakeDataset)

    ts = DatasetSeries(pattern, my_unsupported_kwarg=None)

    assert_raises(TypeError, ts.__getitem__, 0)


def test_get_by_key(FakeDataset, fake_datasets):
    file_list, _pfile_list, _sfile_list, pattern = fake_datasets
    ts = DatasetSeries(pattern)

    Ntot = len(file_list)

    assert ts[0] == ts.get_by_key("current_time", -1)
    assert ts[0] == ts.get_by_key("current_time", 0)
    assert ts[1] == ts.get_by_key("current_time", 0.8)
    assert ts[1] == ts.get_by_key("current_time", 1.2)
    assert ts[Ntot - 1] == ts.get_by_key("current_time", Ntot - 1)
    assert ts[Ntot - 1] == ts.get_by_key("current_time", Ntot)

    assert ts[1] == ts.get_by_key("current_opposite_time", -1.2)
    assert ts[1] == ts.get_by_key("current_opposite_time", -1)
    assert ts[1] == ts.get_by_key("current_opposite_time", -0.6)

import tempfile
from pathlib import Path

import pytest

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
        with pytest.raises(FileNotFoundError):
            DatasetSeries._get_filenames_from_glob_pattern(pattern)


@pytest.fixture
def FakeDataset():
    class _FakeDataset(Dataset):
        """A minimal loadable fake dataset subclass"""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        @classmethod
        def _is_valid(cls, *args, **kwargs):
            return True

        def _parse_parameter_file(self):
            return

        def _set_code_unit_attributes(self):
            return

        def set_code_units(self):
            i = int(Path(self.filename).name.split("_")[-1])
            self.current_time = i
            self.current_redshift = 1 / (i + 1)
            return

        def _hash(self):
            return

        def _setup_classes(self):
            return

    try:
        yield _FakeDataset
    finally:
        output_type_registry.pop("_FakeDataset")


@pytest.fixture
def fake_datasets():
    file_list = [f"fake_data_file_{i:04d}" for i in range(10)]
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
    with pytest.raises(FileNotFoundError):
        DatasetSeries(str(file_list))

    # finally, check that ts[0] fails to actually load
    with pytest.raises(YTUnidentifiedDataType):
        ts[0]


def test_init_fake_dataseries2(FakeDataset, fake_datasets):
    _file_list, _pfile_list, _sfile_list, pattern = fake_datasets
    ds = DatasetSeries(pattern)[0]
    assert isinstance(ds, FakeDataset)

    ts = DatasetSeries(pattern, my_unsupported_kwarg=None)

    with pytest.raises(TypeError):
        ts[0]


def test_get_by_key(FakeDataset, fake_datasets):
    _file_list, _pfile_list, sfile_list, pattern = fake_datasets
    ts = DatasetSeries(pattern)

    Ntot = len(sfile_list)

    t = ts[0].quan(1, "code_time")

    assert sfile_list[0] == ts.get_by_time(-t).filename
    assert sfile_list[0] == ts.get_by_time(t - t).filename
    assert sfile_list[1] == ts.get_by_time((0.8, "code_time")).filename
    assert sfile_list[1] == ts.get_by_time((1.2, "code_time")).filename
    assert sfile_list[Ntot - 1] == ts.get_by_time(t * (Ntot - 1)).filename
    assert sfile_list[Ntot - 1] == ts.get_by_time(t * Ntot).filename

    with pytest.raises(ValueError):
        ts.get_by_time(-2 * t, tolerance=0.1 * t)
    with pytest.raises(ValueError):
        ts.get_by_time(1000 * t, tolerance=0.1 * t)

    assert sfile_list[1] == ts.get_by_redshift(1 / 2.2).filename
    assert sfile_list[1] == ts.get_by_redshift(1 / 2).filename
    assert sfile_list[1] == ts.get_by_redshift(1 / 1.6).filename

    with pytest.raises(ValueError):
        ts.get_by_redshift(1000, tolerance=0.1)

    zmid = (ts[0].current_redshift + ts[1].current_redshift) / 2

    assert sfile_list[1] == ts.get_by_redshift(zmid, prefer="smaller").filename
    assert sfile_list[0] == ts.get_by_redshift(zmid, prefer="larger").filename

import tempfile
from pathlib import Path

from yt.data_objects.time_series import DatasetSeries
from yt.testing import assert_raises
from yt.utilities.exceptions import YTUnidentifiedDataType


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


def test_init_fake_dataseries():

    file_list = [f"fake_data_file_{str(i).zfill(4)}" for i in range(10)]
    with tempfile.TemporaryDirectory() as tmpdir:
        pfile_list = [Path(tmpdir) / file for file in file_list]
        sfile_list = [str(file) for file in pfile_list]
        for file in pfile_list:
            file.touch()
        pattern = Path(tmpdir) / "fake_data_file_*"

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

import os
import tempfile
from pathlib import Path

from yt.data_objects.time_series import get_filenames_from_glob_pattern
from yt.testing import assert_raises


def test_pattern_expansion():
    file_list = [f"fake_data_file_{str(i).zfill(4)}" for i in range(10)]

    with tempfile.TemporaryDirectory() as tmpdir:
        for file in file_list:
            (Path(tmpdir) / file).touch()

        pattern = os.path.join(tmpdir, "fake_data_file_*")
        found = get_filenames_from_glob_pattern(pattern)
        assert found == [os.path.join(tmpdir, file) for file in file_list]

        found2 = get_filenames_from_glob_pattern(Path(pattern))
        assert found2 == found


def test_no_match_pattern():
    with tempfile.TemporaryDirectory() as tmpdir:
        pattern = os.path.join(tmpdir, "fake_data_file_*")
        assert_raises(OSError, get_filenames_from_glob_pattern, pattern)

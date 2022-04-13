import re
import tempfile
from collections import namedtuple
from itertools import chain
from pathlib import Path

import pytest

from yt.frontends.ramses.data_structures import RAMSESFileSanitizer

PathTuple = namedtuple(
    "PathTuple", ("output_dir", "group_dir_name", "info_file", "paths_to_try")
)


def generate_paths(create):
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "output_00123"
        output_dir.mkdir()
        amr_file = output_dir / "amr_00123.out00001"
        if create:
            amr_file.touch()
        info_file = output_dir / "info_00123.txt"
        if create:
            info_file.touch()

        # Test with regular structure
        output_dir2 = Path(tmpdir) / "output_00124"
        output_dir2.mkdir()

        group_dir2 = output_dir2 / "group_00001"
        group_dir2.mkdir()
        info_file2 = group_dir2 / "info_00124.txt"
        if create:
            info_file2.touch()
        amr_file2 = group_dir2 / "amr_00124.out00001"
        if create:
            amr_file2.touch()

        yield (
            PathTuple(
                output_dir=output_dir,
                group_dir_name=None,
                info_file=info_file,
                paths_to_try=(output_dir, info_file, amr_file),
            ),
            PathTuple(
                output_dir=output_dir2,
                group_dir_name=group_dir2.name,
                info_file=info_file2,
                paths_to_try=(output_dir2, info_file2, group_dir2, amr_file2),
            ),
        )


@pytest.fixture
def valid_path_tuples():
    yield from generate_paths(create=True)


@pytest.fixture
def invalid_path_tuples():
    yield from generate_paths(create=False)


def test_valid_sanitizing(valid_path_tuples):
    for answer in valid_path_tuples:
        for path in answer.paths_to_try:
            sanitizer = RAMSESFileSanitizer(path)
            sanitizer.validate()

            assert sanitizer.root_folder == answer.output_dir
            assert sanitizer.group_name == answer.group_dir_name
            assert sanitizer.info_fname == answer.info_file


def test_invalid_sanitizing(valid_path_tuples, invalid_path_tuples):
    for path in chain(*(pt.paths_to_try for pt in invalid_path_tuples)):
        sanitizer = RAMSESFileSanitizer(path)

        if path.exists():
            expected_error = ValueError
            expected_error_message = (
                "Could not determine output directory from '.*'\n"
                "Expected a directory name of form .* "
                "containing an info_\\*.txt file and amr_\\* files."
            )
        else:
            expected_error = FileNotFoundError
            expected_error_message = re.escape(
                f"No such file or directory '{str(path)}'"
            )

        with pytest.raises(expected_error, match=expected_error_message):
            sanitizer.validate()

    for path in chain(*(pt.paths_to_try for pt in valid_path_tuples)):
        expected_error_message = re.escape(
            f"No such file or directory '{str(path/'does_not_exist.txt')}'"
        )
        sanitizer = RAMSESFileSanitizer(path / "does_not_exist.txt")
        with pytest.raises(FileNotFoundError, match=expected_error_message):
            sanitizer.validate()

    expected_error_message = "No such file or directory '.*'"
    sanitizer = RAMSESFileSanitizer(Path("this") / "does" / "not" / "exist")
    with pytest.raises(FileNotFoundError, match=expected_error_message):
        sanitizer.validate()

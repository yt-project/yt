import numpy as np
import pytest

from yt.utilities.cython_fortran_utils import FortranFile


def test_raise_error_when_file_does_not_exist():
    with pytest.raises(FileNotFoundError):
        FortranFile("/this/file/does/not/exist")


def test_read(tmp_path):
    dummy_file = tmp_path / "test.bin"
    # Write a Fortran-formatted file containing one record with 4 doubles
    dummy_file.write_bytes(
        b" \x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x08@\x00\x00\x00\x00\x00\x00\x10@ \x00\x00\x00"
    )
    with FortranFile(str(dummy_file)) as f:
        np.testing.assert_equal(
            f.read_vector("d"),
            [1.0, 2.0, 3.0, 4.0],
        )

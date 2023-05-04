import os

import numpy as np
import pytest

from yt.utilities.cython_fortran_utils import FortranFile


def test_raise_error_when_file_does_not_exist():
    with pytest.raises(FileNotFoundError):
        FortranFile("/this/file/does/not/exist")


def test_read(tmpdir):
    dummy_file = os.path.join(tmpdir, "test.bin")
    # Write a Fortran-formatted file containing one record with 4 doubles
    with open(dummy_file, "bw") as f:
        f.write(
            b" \x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x08@\x00\x00\x00\x00\x00\x00\x10@ \x00\x00\x00"
        )

    with FortranFile(dummy_file) as f:
        np.testing.assert_equal(
            f.read_vector("d"),
            [1.0, 2.0, 3.0, 4.0],
        )

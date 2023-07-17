import struct

import numpy as np
import pytest

from yt.utilities.cython_fortran_utils import FortranFile


def test_raise_error_when_file_does_not_exist():
    with pytest.raises(FileNotFoundError):
        FortranFile("/this/file/does/not/exist")


def test_read(tmp_path):
    dummy_file = tmp_path / "test.bin"
    # Write a Fortran-formatted file containing one record with 4 doubles
    # The format is a 32bit integer with value 4*sizeof(double)=32
    # followed by 4 doubles and another 32bit integer with value 32
    # Note that there is no memory alignment, hence the "=" below
    buff = struct.pack("=i 4d i", 32, 1.0, 2.0, 3.0, 4.0, 32)
    dummy_file.write_bytes(buff)
    with FortranFile(str(dummy_file)) as f:
        np.testing.assert_equal(
            f.read_vector("d"),
            [1.0, 2.0, 3.0, 4.0],
        )

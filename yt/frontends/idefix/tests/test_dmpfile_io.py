import os

import numpy as np
from numpy.testing import assert_equal

from yt.frontends.idefix.dmpfile_io import read_idefix_dmpfile
from yt.testing import requires_file

test_file = os.path.join(os.environ["IDEFIX_DIR"], "test", "HD", "KHI", "dump.0001.dmp")


@requires_file(test_file)
def test_read_dmp():
    with open(test_file, mode="rb") as fh:
        fprops, fdata = read_idefix_dmpfile(fh)

    expected_fields = [
        "x1",
        "x2",
        "x3",
        "Vc-RHO",
        "Vc-VX1",
        "Vc-VX2",
        "time",
        "dt",
        "vtkFileNumber",
        "vtktnext",
        "dumpFileNumber",
        "dumptnext",
    ]
    detected_fields = list(fprops.keys())
    assert detected_fields == expected_fields

    expected_shape = [1024, 256, 1]
    detected_shape = np.concatenate([fprops[k][-1] for k in ("x1", "x2", "x3")])
    assert_equal(detected_shape, expected_shape)
    for field_name, data in fprops.items():
        if not field_name.startswith("V"):
            continue
        _dtype, _ndim, dim = data
        assert_equal(dim, expected_shape)
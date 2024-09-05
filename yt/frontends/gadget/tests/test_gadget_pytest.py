import numpy as np

import yt
from yt.testing import requires_file, requires_module
from yt.utilities.on_demand_imports import _h5py as h5py


@requires_file("snapshot_033/snap_033.0.hdf5")
@requires_module("h5py")
def test_gadget_header_array_reduction(tmp_path):
    # first get a real header
    ds = yt.load("snapshot_033/snap_033.0.hdf5")
    hvals = ds._get_hvals()
    hvals_orig = hvals.copy()
    # wrap some of the scalar values in nested arrays
    hvals["Redshift"] = np.array([hvals["Redshift"]])
    hvals["Omega0"] = np.array([[hvals["Omega0"]]])

    # drop those header values into a fake header-only file
    tmp_snpshot_dir = tmp_path / "snapshot_033"
    tmp_snpshot_dir.mkdir()
    tmp_header_only_file = str(tmp_snpshot_dir / "fake_gadget_header.hdf5")
    with h5py.File(tmp_header_only_file, mode="w") as f:
        headergrp = f.create_group("Header")
        for field, val in hvals.items():
            headergrp.attrs[field] = val

    # trick the dataset into using the header file and make sure the
    # arrays are reduced
    ds._input_filename = tmp_header_only_file
    hvals = ds._get_hvals()
    for attr in ("Redshift", "Omega0"):
        assert hvals[attr] == hvals_orig[attr]
        assert isinstance(hvals[attr], np.ndarray) is False

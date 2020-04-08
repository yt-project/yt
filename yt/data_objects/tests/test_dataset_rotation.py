# writeme

import yt
from yt.testing import fake_amr_ds
import numpy as np


def count_invalid_values(arr):
    return len(arr[~np.isfinite(arr)])

def get_buffer(ds):
    slc = ds.r[:,0.5,:] # zslice
    buff = ds.coordinates.pixelize(dimension=1, data_source=slc,
            field=('stream', 'Density'), bounds=[-1., 1., -1., 1.], size=(800, 800))
    return buff


def test_dont_corrupt_data():
    ds = fake_amr_ds(geometry="cylindrical")

    ndiv = 11 # arbitrary, not a power of two on purpose

    buff = get_buffer(ds)
    ref_invalid_count = count_invalid_values(buff)
    assert ref_invalid_count < 800*800

    for i in range(ndiv+1):
        ds.rotate_frame_by(-np.pi/ndiv)
        slc = ds.r[:,0.5,:] # zslice
        assert count_invalid_values(slc["Density"]) == 0

        buff = get_buffer(ds)
        new_invalid_count = count_invalid_values(buff)
        assert new_invalid_count == ref_invalid_count
    
        slc.to_pw()

def test_dont_corrupt_data_SlicePlots():
    ds = fake_amr_ds(geometry="cylindrical")

    ndiv = 11 # arbitrary, not a power of two on purpose

    buff = get_buffer(ds)
    ref_invalid_count = count_invalid_values(buff)
    assert ref_invalid_count < 800*800

    for i in range(ndiv+1):
        ds.rotate_frame_by(-np.pi/ndiv)
        slc = ds.r[:,0.5,:] # zslice
        assert count_invalid_values(slc["Density"]) == 0

        buff = get_buffer(ds)
        new_invalid_count = count_invalid_values(buff)
        assert new_invalid_count == ref_invalid_count

        yt.SlicePlot(ds, "z", "Density") # <- THIS is the line that actually corrupts data, so it seems
# writeme

from itertools import product
import yt
from yt.testing import fake_amr_ds, assert_allclose
import numpy as np

import pytest

def count_invalid_values(arr):
    return len(arr[~np.isfinite(arr)])

def get_buffer(ds):
    slc = ds.r[:,0.5,:] # zslice
    buff = ds.coordinates.pixelize(dimension=1, data_source=slc,
            field=('stream', 'Density'), bounds=[-1., 1., -1., 1.], size=(800, 800))
    return buff

def plot_topw(ds):
    slc = ds.r[:,0.5,:]
    return slc.to_pw('Density')

def plot_SlicePlot(ds):
    return yt.SlicePlot(ds, 'z', 'Density')

def inf2zero(arr):
    arr[arr == np.inf] = 0

methods = [lambda x: None, plot_topw, plot_SlicePlot]

@pytest.mark.parametrize('plot_method,touch', list(product(methods, [True, False])))
def test_dont_corrupt_data(plot_method, touch):
    ds = fake_amr_ds(geometry="cylindrical")

    ndiv = 3 # arbitrary, not a power of two on purpose

    if touch:
        # devnote : I do not understand why, but creating a buffer NOW affects the end result
        # later when convert a rotated slice to a plot,
        # however this only works in combination with plotting_method=plot_topw
        buff0 = get_buffer(ds)

        # invalid values (nan, inf) are used to fill the buffer in the absence of data,
        # typically in the corners of a plot window for cylindrical data seen from the z axis
        ref_invalid_count = count_invalid_values(buff0)
        assert ref_invalid_count < 800*800

        inf2zero(buff0)
        # the minimal code with identical effect simply consist in touching (hence caching ?)
        # the pixel-x coord in the slice object
        """
        slc = ds.r[:,0.5,:] # zslice
        slc["px"]
        """

    for i in range(ndiv+1):
        ds.rotate_frame_by(-np.pi/ndiv, safe_=False)
        slc = ds.r[:,0.5,:] # zslice
        assert count_invalid_values(slc["Density"]) == 0
        
        if touch:
            buff1 = get_buffer(ds)
            new_invalid_count = count_invalid_values(buff1)
            # although the exact number of filler pixels may slighlty vary, a variation greater than 0.1%
            # is very suspicious in this test case
            assert_allclose(new_invalid_count, ref_invalid_count, rtol=1e-3)
            inf2zero(buff1)
            assert_allclose(buff1.mean(), buff0.mean(), rtol=1e-3)
        

        plot = plot_method(ds)
        if plot is not None:
            buff3 = plot.frb.data['Density']
            assert count_invalid_values(buff3) < buff3.size/2
        
        
        if touch:
            buff2 = get_buffer(ds)
            new_invalid_count = count_invalid_values(buff2)
            assert_allclose(new_invalid_count, ref_invalid_count, rtol=1e-3)
            inf2zero(buff2)
            assert_allclose(buff1, buff2, rtol=1e-14)
        if plot is not None:
            buff3 = plot.frb.data['Density']
            inf2zero(buff3)

            if touch:
                assert_allclose(buff1, buff3, rtol=1e-14)
            #plot.save(f'/tmp/{plot_method.__name__}_{touch}_{i}')
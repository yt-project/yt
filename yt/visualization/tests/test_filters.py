#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
#-----------------------------------------------------------------------------

"""
Tests for frb filters

"""

from yt.testing import fake_amr_ds, requires_module
import yt
import numpy as np

@requires_module("scipy")
def test_white_noise_filter():
    ds = fake_amr_ds(fields=("density",))
    p = ds.proj("density", "z")
    frb = p.to_frb((1, 'unitary'), 64)
    frb.apply_white_noise()
    frb.apply_white_noise(1e-3)
    frb["density"]


@requires_module("scipy")
def test_gauss_beam_filter():
    ds = fake_amr_ds(fields=("density",))
    p = ds.proj("density", "z")
    frb = p.to_frb((1, 'unitary'), 64)
    frb.apply_gauss_beam(sigma=1.0)
    frb["density"]


def test_custom_filter():
    ds = fake_amr_ds(fields=("density",))
    p = ds.proj("density", "z")
    frb = p.to_frb((1, 'unitary'), 64)
    data1 = frb['density']
    frb.apply_custom_filter(lambda d: d*2)
    data2 = frb['density']

    assert np.allclose(data1.value * 2, data2.value)


@requires_module("scipy")
def test_filter_wiring():
    from scipy.ndimage import gaussian_filter
    ds = fake_amr_ds(fields=("density", ))
    p = yt.SlicePlot(ds, 'x', 'density')

    frb1 = p.frb
    data_orig = frb1['density'].value

    sigma = 2
    p.frb.apply_gauss_beam(sigma)
    frb2 = p.frb
    data_gauss = frb2['density'].value

    p.frb.apply_white_noise()
    frb3 = p.frb
    data_white = frb3['density'].value

    # We check the frb have changed
    assert frb1 != frb2
    assert frb1 != frb3
    assert frb2 != frb3

    # We check the resulting image are different each time
    assert np.allclose(data_orig, data_gauss) is False
    assert np.allclose(data_orig, data_white) is False
    assert np.allclose(data_gauss, data_white) is False

    # Check the gaussian filtering should is ok
    assert np.allclose(gaussian_filter(data_orig, sigma),
                       data_gauss)

"""
Tests for frb filters

"""

from yt.testing import fake_amr_ds, requires_module
import yt
import numpy as np

@requires_module("scipy")
def test_white_noise_filter():
    ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
    p = ds.proj(("gas", "density"), "z")
    frb = p.to_frb((1, "unitary"), 64)
    frb.apply_white_noise()
    frb.apply_white_noise(1e-3)
    frb[("gas", "density")]


@requires_module("scipy")
def test_gauss_beam_filter():
    ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
    p = ds.proj(("gas", "density"), "z")
    frb = p.to_frb((1, "unitary"), 64)
    frb.apply_gauss_beam(sigma=1.0)
    frb[("gas", "density")]


@requires_module("scipy")
def test_filter_wiring():
    from scipy.ndimage import gaussian_filter
    ds = fake_amr_ds(fields=("density", ))
    p = yt.SlicePlot(ds, 'x', 'density')
    data_before = p.frb['density'].value
    sigma = 2
    p.frb.apply_gauss_beam(sigma)
    data_after = p.frb['density'].value

    # Check that the frb yields the same result
    assert np.allclose(gaussian_filter(data_before, sigma),
                       data_after)

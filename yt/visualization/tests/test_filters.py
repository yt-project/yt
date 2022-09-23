"""
Tests for frb filters

"""

import numpy as np

import yt
from yt.testing import fake_amr_ds, requires_module


@requires_module("scipy")
def test_white_noise_filter():
    ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
    p = ds.proj(("gas", "density"), "z")
    frb = p.to_frb((1, "unitary"), 64)
    frb.apply_white_noise()
    frb.apply_white_noise(1e-3)
    frb.render(("gas", "density"))


@requires_module("scipy")
def test_gauss_beam_filter():
    ds = fake_amr_ds(fields=("density",), units=("g/cm**3",))
    p = ds.proj(("gas", "density"), "z")
    frb = p.to_frb((1, "unitary"), 64)
    frb.apply_gauss_beam(sigma=1.0)
    frb.render(("gas", "density"))


@requires_module("scipy")
def test_filter_wiring():

    ds = fake_amr_ds(fields=[("gas", "density")], units=["g/cm**3"])
    p = yt.SlicePlot(ds, "x", "density")

    # Note: frb is a FixedResolutionBuffer object
    frb1 = p.frb
    data_orig = frb1["density"].value

    sigma = 2
    nbeam = 30
    p.frb.apply_gauss_beam(nbeam=nbeam, sigma=sigma)
    frb2 = p.frb
    data_gauss = frb2["density"].value

    p.frb.apply_white_noise()
    frb3 = p.frb
    data_white = frb3["density"].value

    # We check the frb objects are different
    assert frb1 is not frb2
    assert frb1 is not frb3
    assert frb2 is not frb3

    # We check the resulting image are different each time
    assert not np.allclose(data_orig, data_gauss)
    assert not np.allclose(data_orig, data_white)
    assert not np.allclose(data_gauss, data_white)

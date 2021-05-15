"""
Tests for frb filters

"""

from yt.testing import fake_amr_ds, requires_module


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
    frb.apply_gauss_beam(nbeam=15, sigma=1.0)
    frb[("gas", "density")]

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
    frb.apply_gauss_beam(nbeam=15, sigma=1.0)
    frb["density"]

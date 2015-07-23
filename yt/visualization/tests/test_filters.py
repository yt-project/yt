#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
#-----------------------------------------------------------------------------

"""
Tests for frb filters

"""

from yt.testing import fake_amr_ds


class TestFilters():

    @classmethod
    def setup_class(cls):
        ds = fake_amr_ds(fields=("density",))
        p = ds.proj("density", "z")
        cls.frb = p.to_frb((1, 'unitary'), 64)

    def teardown(self):
        try:
            del self.frb["density"]
        except KeyError:
            pass

    def test_white_noise_filter(self):
        self.frb.apply_white_noise()
        self.frb.apply_white_noise(1e-3)
        self.frb["density"]

    def test_gauss_beam_filter(self):
        self.frb.apply_gauss_beam(nbeam=15, sigma=1.0)
        self.frb["density"]

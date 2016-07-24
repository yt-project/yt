"""
GAMER frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    data_dir_load
from yt.frontends.gamer.api import GAMERDataset



jet         = "InteractingJets/jet_000002"
_fields_jet = ("temperature", "density", "velocity_magnitude")
jet_units   = {"length_unit":(1.0,"kpc"),
               "time_unit"  :(3.08567758096e+13,"s"),
               "mass_unit"  :(1.4690033e+36,"g")}

@requires_ds(jet, big_data=True)
def test_jet():
    ds = data_dir_load(jet, kwargs={"units_override":jet_units})
    yield assert_equal, str(ds), "jet_000002"
    for test in small_patch_amr(ds, _fields_jet):
        test_jet.__name__ = test.description
        yield test


psiDM         = "WaveDarkMatter/psiDM_000020"
_fields_psiDM = ("Dens", "Real", "Imag")

@requires_ds(psiDM, big_data=True)
def test_psiDM():
    ds = data_dir_load(psiDM)
    yield assert_equal, str(ds), "psiDM_000020"
    for test in small_patch_amr(ds, _fields_psiDM):
        test_psiDM.__name__ = test.description
        yield test


@requires_file(psiDM)
def test_GAMERDataset():
    assert isinstance(data_dir_load(psiDM), GAMERDataset)


@requires_file(jet)
def test_units_override():
    for test in units_override_check(jet):
        yield test

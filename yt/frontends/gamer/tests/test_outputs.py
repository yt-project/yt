from yt.frontends.gamer.api import GAMERDataset
from yt.testing import assert_equal, requires_file, units_override_check
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

jet = "InteractingJets/jet_000002"
_fields_jet = (
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
)
jet_units = {
    "length_unit": (1.0, "kpc"),
    "time_unit": (3.08567758096e13, "s"),
    "mass_unit": (1.4690033e36, "g"),
}


@requires_ds(jet, big_data=True)
def test_jet():
    ds = data_dir_load(jet, kwargs={"units_override": jet_units})
    assert_equal(str(ds), "jet_000002")
    for test in small_patch_amr(ds, _fields_jet):
        test_jet.__name__ = test.description
        yield test


psiDM = "WaveDarkMatter/psiDM_000020"
_fields_psiDM = ("Dens", "Real", "Imag")


@requires_ds(psiDM, big_data=True)
def test_psiDM():
    ds = data_dir_load(psiDM)
    assert_equal(str(ds), "psiDM_000020")
    for test in small_patch_amr(ds, _fields_psiDM):
        test_psiDM.__name__ = test.description
        yield test


plummer = "Plummer/plummer_000000"
_fields_plummer = (("gamer", "ParDens"), ("deposit", "io_cic"))


@requires_ds(plummer, big_data=True)
def test_plummer():
    ds = data_dir_load(plummer)
    assert_equal(str(ds), "plummer_000000")
    for test in small_patch_amr(ds, _fields_plummer):
        test_plummer.__name__ = test.description
        yield test


mhdvortex = "MHDOrszagTangVortex/Data_000018"
_fields_mhdvortex = (
    ("gamer", "CCMagX"),
    ("gamer", "CCMagY"),
    ("gas", "magnetic_energy"),
)


@requires_ds(mhdvortex, big_data=True)
def test_mhdvortex():
    ds = data_dir_load(mhdvortex)
    assert_equal(str(ds), "Data_000018")
    for test in small_patch_amr(ds, _fields_mhdvortex):
        test_mhdvortex.__name__ = test.description
        yield test


@requires_file(psiDM)
def test_GAMERDataset():
    assert isinstance(data_dir_load(psiDM), GAMERDataset)


@requires_file(jet)
def test_units_override():
    units_override_check(jet)


jiw = "JetICMWall/Data_000060"
_fields_jiw = (
    ("gas", "four_velocity_magnitude"),
    ("gas", "density"),
    ("gas", "gamma"),
    ("gas", "temperature"),
)


@requires_ds(jiw, big_data=True)
def test_jiw():
    ds = data_dir_load(jiw)
    assert_equal(str(ds), "Data_000060")
    for test in small_patch_amr(ds, _fields_jiw):
        test_jiw.__name__ = test.description
        yield test

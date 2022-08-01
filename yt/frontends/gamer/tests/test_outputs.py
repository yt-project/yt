from yt.frontends.gamer.api import GAMERDataset
from yt.testing import (
    assert_array_almost_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
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
    ("gas", "magnetic_energy_density"),
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


@requires_ds(jiw, big_data=True)
def test_stress_energy():
    axes = "txyz"
    ds = data_dir_load(jiw)
    center = ds.arr([3.0, 10.0, 10.0], "kpc")
    sp = ds.sphere(center, (1.0, "kpc"))
    c2 = ds.units.clight**2
    inv_c2 = 1.0 / c2
    rho = sp["gas", "density"]
    p = sp["gas", "pressure"]
    e = sp["gas", "thermal_energy_density"]
    h = rho + (e + p) * inv_c2
    for mu in range(4):
        for nu in range(4):
            # matrix is symmetric so only do the upper-right part
            if nu >= mu:
                Umu = sp[f"four_velocity_{axes[mu]}"]
                Unu = sp[f"four_velocity_{axes[nu]}"]
                Tmunu = h * Umu * Unu
                if mu != nu:
                    assert_array_almost_equal(sp[f"T{mu}{nu}"], sp[f"T{nu}{mu}"])
                else:
                    Tmunu += p
                assert_array_almost_equal(sp[f"T{mu}{nu}"], Tmunu)

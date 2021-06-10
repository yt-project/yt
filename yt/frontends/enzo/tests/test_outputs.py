import numpy as np

from yt.frontends.enzo.api import EnzoDataset
from yt.frontends.enzo.fields import NODAL_FLAGS
from yt.testing import (
    assert_allclose_units,
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    big_patch_amr,
    data_dir_load,
    requires_ds,
    small_patch_amr,
)
from yt.visualization.plot_window import SlicePlot

_fields = (
    ("gas", "temperature"),
    ("gas", "density"),
    ("gas", "velocity_magnitude"),
    ("gas", "velocity_divergence"),
)

two_sphere_test = "ActiveParticleTwoSphere/DD0011/DD0011"
active_particle_cosmology = "ActiveParticleCosmology/DD0046/DD0046"
ecp = "enzo_cosmology_plus/DD0046/DD0046"
g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
toro1d = "ToroShockTube/DD0001/data0001"
kh2d = "EnzoKelvinHelmholtz/DD0011/DD0011"
mhdctot = "MHDCTOrszagTang/DD0004/data0004"
dnz = "DeeplyNestedZoom/DD0025/data0025"
p3mini = "PopIII_mini/DD0034/DD0034"


def color_conservation(ds):
    species_names = ds.field_info.species_names
    dd = ds.all_data()
    dens_yt = dd[("gas", "density")].copy()
    # Enumerate our species here
    for s in sorted(species_names):
        if s == "El":
            continue
        dens_yt -= dd[f"{s}_density"]
    dens_yt -= dd[("enzo", "Metal_Density")]
    delta_yt = np.abs(dens_yt / dd[("gas", "density")])
    # Now we compare color conservation to Enzo's color conservation
    dd = ds.all_data()
    dens_enzo = dd[("enzo", "Density")].copy()
    for f in sorted(ds.field_list):
        ff = f[1]
        if not ff.endswith("_Density"):
            continue
        start_strings = [
            "Electron_",
            "SFR_",
            "Forming_Stellar_",
            "Dark_Matter",
            "Star_Particle_",
        ]
        if any([ff.startswith(ss) for ss in start_strings]):
            continue
        dens_enzo -= dd[f]
    delta_enzo = np.abs(dens_enzo / dd[("enzo", "Density")])
    np.testing.assert_almost_equal(delta_yt, delta_enzo)


def check_color_conservation(ds):
    species_names = ds.field_info.species_names
    dd = ds.all_data()
    dens_yt = dd[("gas", "density")].copy()
    # Enumerate our species here
    for s in sorted(species_names):
        if s == "El":
            continue
        dens_yt -= dd[f"{s}_density"]
    dens_yt -= dd[("enzo", "Metal_Density")]
    delta_yt = np.abs(dens_yt / dd[("gas", "density")])

    # Now we compare color conservation to Enzo's color conservation
    dd = ds.all_data()
    dens_enzo = dd[("enzo", "Density")].copy()
    for f in sorted(ds.field_list):
        ff = f[1]
        if not ff.endswith("_Density"):
            continue
        start_strings = [
            "Electron_",
            "SFR_",
            "Forming_Stellar_",
            "Dark_Matter",
            "Star_Particle_",
        ]
        if any([ff.startswith(ss) for ss in start_strings]):
            continue
        dens_enzo -= dd[f]
    delta_enzo = np.abs(dens_enzo / dd[("enzo", "Density")])
    return assert_almost_equal, delta_yt, delta_enzo


m7 = "DD0010/moving7_0010"


@requires_ds(m7)
def test_moving7():
    ds = data_dir_load(m7)
    assert_equal(str(ds), "moving7_0010")
    for test in small_patch_amr(m7, _fields):
        test_moving7.__name__ = test.description
        yield test


@requires_ds(g30, big_data=True)
def test_galaxy0030():
    ds = data_dir_load(g30)
    yield check_color_conservation(ds)
    assert_equal(str(ds), "galaxy0030")
    for test in big_patch_amr(ds, _fields):
        test_galaxy0030.__name__ = test.description
        yield test
    assert_equal(ds.particle_type_counts, {"io": 1124453})


@requires_ds(toro1d)
def test_toro1d():
    ds = data_dir_load(toro1d)
    assert_equal(str(ds), "data0001")
    for test in small_patch_amr(ds, ds.field_list):
        test_toro1d.__name__ = test.description
        yield test


@requires_ds(kh2d)
def test_kh2d():
    ds = data_dir_load(kh2d)
    assert_equal(str(ds), "DD0011")
    for test in small_patch_amr(ds, ds.field_list):
        test_kh2d.__name__ = test.description
        yield test


@requires_ds(ecp, big_data=True)
def test_ecp():
    ds = data_dir_load(ecp)
    # Now we test our species fields
    yield check_color_conservation(ds)


@requires_file(enzotiny)
def test_units_override():
    units_override_check(enzotiny)


@requires_ds(ecp, big_data=True)
def test_nuclei_density_fields():
    ds = data_dir_load(ecp)
    ad = ds.all_data()
    assert_array_equal(
        ad[("gas", "H_nuclei_density")],
        (ad[("gas", "H_p0_number_density")] + ad[("gas", "H_p1_number_density")]),
    )
    assert_array_equal(
        ad[("gas", "He_nuclei_density")],
        (
            ad[("gas", "He_p0_number_density")]
            + ad[("gas", "He_p1_number_density")]
            + ad[("gas", "He_p2_number_density")]
        ),
    )


@requires_file(enzotiny)
def test_EnzoDataset():
    assert isinstance(data_dir_load(enzotiny), EnzoDataset)


@requires_file(two_sphere_test)
@requires_file(active_particle_cosmology)
def test_active_particle_datasets():
    two_sph = data_dir_load(two_sphere_test)
    assert "AccretingParticle" in two_sph.particle_types_raw
    assert "io" not in two_sph.particle_types_raw
    assert "all" in two_sph.particle_types
    assert "nbody" in two_sph.particle_types
    assert_equal(len(two_sph.particle_unions), 2)
    pfields = [
        "GridID",
        "creation_time",
        "dynamical_time",
        "identifier",
        "level",
        "metallicity",
        "particle_mass",
    ]
    pfields += [f"particle_position_{d}" for d in "xyz"]
    pfields += [f"particle_velocity_{d}" for d in "xyz"]

    acc_part_fields = [("AccretingParticle", pf) for pf in ["AccretionRate"] + pfields]

    real_acc_part_fields = sorted(
        f for f in two_sph.field_list if f[0] == "AccretingParticle"
    )
    assert_equal(acc_part_fields, real_acc_part_fields)

    apcos = data_dir_load(active_particle_cosmology)
    assert_equal(["CenOstriker", "DarkMatter"], apcos.particle_types_raw)
    assert "all" in apcos.particle_unions
    assert "nbody" in apcos.particle_unions

    apcos_fields = [("CenOstriker", pf) for pf in pfields]

    real_apcos_fields = sorted(f for f in apcos.field_list if f[0] == "CenOstriker")

    assert_equal(apcos_fields, real_apcos_fields)

    assert_equal(
        apcos.particle_type_counts, {"CenOstriker": 899755, "DarkMatter": 32768}
    )


@requires_file(mhdctot)
def test_face_centered_mhdct_fields():
    ds = data_dir_load(mhdctot)

    ad = ds.all_data()
    grid = ds.index.grids[0]

    for field, flag in NODAL_FLAGS.items():
        dims = ds.domain_dimensions
        assert_equal(ad[field].shape, (dims.prod(), 2 * sum(flag)))
        assert_equal(grid[field].shape, tuple(dims) + (2 * sum(flag),))

    # Average of face-centered fields should be the same as cell-centered field
    assert (ad[("enzo", "BxF")].sum(axis=-1) / 2 == ad[("enzo", "Bx")]).all()
    assert (ad[("enzo", "ByF")].sum(axis=-1) / 2 == ad[("enzo", "By")]).all()
    assert (ad[("enzo", "BzF")].sum(axis=-1) / 2 == ad[("enzo", "Bz")]).all()


@requires_file(dnz)
def test_deeply_nested_zoom():
    ds = data_dir_load(dnz)

    # carefully chosen to just barely miss a grid in the middle of the image
    center = [0.4915073260199302, 0.5052605316800006, 0.4905805557500548]

    plot = SlicePlot(ds, "z", "density", width=(0.001, "pc"), center=center)

    image = plot.frb[("gas", "density")]

    assert (image > 0).all()

    v, c = ds.find_max(("gas", "density"))

    assert_allclose_units(v, ds.quan(0.005878286377124154, "g/cm**3"))

    c_actual = [0.49150732540021, 0.505260532936791, 0.49058055816398]
    c_actual = ds.arr(c_actual, "code_length")
    assert_allclose_units(c, c_actual)

    assert_equal(max(g[("gas", "density")].max() for g in ds.index.grids), v)


@requires_file(kh2d)
def test_2d_grid_shape():
    # see issue #1601
    # we want to make sure that accessing data on a grid object
    # returns a 3D array with a dummy dimension.
    ds = data_dir_load(kh2d)
    g = ds.index.grids[1]
    assert g[("gas", "density")].shape == (128, 100, 1)


@requires_file(p3mini)
def test_nonzero_omega_radiation():
    """
    Test support for non-zero omega_radiation cosmologies.
    """
    ds = data_dir_load(p3mini)

    assert_equal(ds.omega_radiation, ds.cosmology.omega_radiation)

    tratio = ds.current_time / ds.cosmology.t_from_z(ds.current_redshift)
    assert_almost_equal(
        tratio,
        1,
        4,
        err_msg="Simulation time not consistent with cosmology calculator.",
    )

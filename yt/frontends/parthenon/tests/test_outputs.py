import numpy as np

from yt.frontends.parthenon.api import ParthenonDataset
from yt.loaders import load
from yt.testing import (
    assert_allclose,
    assert_equal,
    assert_true,
    requires_file,
)
from yt.utilities.answer_testing.framework import (
    GenericArrayTest,
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

_fields_parthenon_advection = (
    ("parthenon", "advected_0_0"),
    ("parthenon", "one_minus_advected"),
    ("parthenon", "one_minus_advected_sq"),
    ("parthenon", "one_minus_sqrt_one_minus_advected_sq_12"),
    ("parthenon", "one_minus_sqrt_one_minus_advected_sq_37"),
)

# Simple 2D test (advected spherical blob) with AMR from the main Parthenon test suite
# adjusted so that x1 != x2.
# Ran with `./example/advection/advection-example -i ../tst/regression/test_suites/output_hdf5/parthinput.advection parthenon/mesh/nx1=128 parthenon/mesh/x1min=-1.0 parthenon/mesh/x1max=1.0 Advection/vx=2`
# on changeset e5059ad
parthenon_advection = "parthenon_advection/advection_2d.out0.final.phdf"


@requires_ds(parthenon_advection)
def test_loading_data():
    ds = data_dir_load(parthenon_advection)
    assert_equal(str(ds), "advection_2d.out0.final")
    dd = ds.all_data()
    # test mesh dims
    vol = np.prod(ds.domain_right_edge - ds.domain_left_edge)
    assert_equal(vol, ds.quan(2.0, "code_length**3"))
    assert_allclose(dd.quantities.total_quantity("cell_volume"), vol)
    # test data
    for field in _fields_parthenon_advection:

        def field_func(name):
            return dd[name]

        yield GenericArrayTest(ds, field_func, args=[field])

    # reading data of two fields and compare against each other (data is squared in output)
    ad = ds.all_data()
    assert_allclose(
        ad[("parthenon", "one_minus_advected")] ** 2.0,
        ad[("parthenon", "one_minus_advected_sq")],
    )

    # check if the peak is in the domain center (and at the highest refinement level)
    dist_of_max_from_center = np.linalg.norm(
        ad.quantities.max_location(("parthenon", "advected_0_0"))[1:] - ds.domain_center
    )

    dx_min, dx_max = ad.quantities.extrema(("index", "dx"))
    dy_min, dy_max = ad.quantities.extrema(("index", "dy"))

    assert_true(dist_of_max_from_center < np.min((dx_min, dy_min)))


# 3D magnetized cluster center from downstream Parthenon code AthenaPK (Restart, Conserveds)
athenapk_cluster = "athenapk_cluster/athenapk_cluster.restart.00000.rhdf"

# Keplerian disk in 2D cylindrical from downstream Parthenon code AthenaPK (Data, Primitives)
athenapk_disk = "athenapk_disk/athenapk_disk.prim.00000.phdf"


@requires_file(athenapk_cluster)
def test_AthenaPK_rhdf():
    # Test that a downstream AthenaPK data set can be loaded with this Parthenon
    # frontend
    ds = data_dir_load(athenapk_cluster)
    assert isinstance(ds, ParthenonDataset)

    assert_equal(ds.domain_left_edge.in_units("code_length").v, (-0.15, -0.18, -0.2))
    assert_equal(ds.domain_right_edge.in_units("code_length").v, (0.15, 0.18, 0.2))


@requires_file(athenapk_disk)
def test_AthenaPK_phdf():
    # Test that a downstream AthenaPK data set can be loaded with this Parthenon
    # frontend
    assert isinstance(data_dir_load(athenapk_disk), ParthenonDataset)


_fields_derived = (
    ("gas", "temperature"),
    ("gas", "specific_thermal_energy"),
)

_fields_derived_cluster = (("gas", "magnetic_field_strength"),)


@requires_ds(athenapk_cluster)
def test_cluster():
    ds = data_dir_load(athenapk_cluster)
    assert_equal(str(ds), "athenapk_cluster.restart.00000")
    for test in small_patch_amr(ds, _fields_derived + _fields_derived_cluster):
        test_cluster.__name__ = test.description
        yield test


@requires_ds(athenapk_disk)
@requires_ds(athenapk_cluster)
def test_derived_fields():
    # Check that derived fields like temperature are present in downstream
    # which define them

    # Check temperature and specific thermal energy becomes defined from primitives
    ds = data_dir_load(athenapk_disk)
    dd = ds.all_data()

    for field in _fields_derived:

        def field_func(name):
            return dd[name]

        yield GenericArrayTest(ds, field_func, args=[field])

    # Check hydro, magnetic, and cooling fields defined from conserveds
    ds = data_dir_load(athenapk_cluster)
    dd = ds.all_data()

    for field in _fields_derived + _fields_derived_cluster:

        def field_func(name):
            return dd[name]

        yield GenericArrayTest(ds, field_func, args=[field])


@requires_file(athenapk_cluster)
@requires_file(athenapk_disk)
def test_adiabatic_index():
    # Read adiabiatic index from dataset parameters
    ds = data_dir_load(athenapk_cluster)
    assert_allclose(ds.gamma, 5.0 / 3.0, rtol=1e-12)

    ds = data_dir_load(athenapk_disk)
    assert_allclose(ds.gamma, 4.0 / 3.0, rtol=1e-12)

    # Change adiabatic index from dataset parameters
    ds = load(athenapk_disk, parameters={"gamma": 9.0 / 8.0})
    assert_allclose(ds.gamma, 9.0 / 8.0, rtol=1e-12)


@requires_file(athenapk_cluster)
def test_molecular_mass():
    # Read mu from dataset parameters
    ds = data_dir_load(athenapk_cluster)
    assert_allclose(float(ds.mu), 0.5925925925925926, rtol=1e-12)

    # Change He mass fraction from dataset parameters
    ds = load(athenapk_disk, parameters={"mu": 137})
    assert_equal(ds.mu, 137)


@requires_file(athenapk_cluster)
def test_units():
    # Check units in dataset are loaded correctly
    ds = data_dir_load(athenapk_cluster)
    assert_allclose(float(ds.quan(1, "code_time").in_units("Gyr")), 1, rtol=1e-12)
    assert_allclose(float(ds.quan(1, "code_length").in_units("Mpc")), 1, rtol=1e-12)
    assert_allclose(float(ds.quan(1, "code_mass").in_units("msun")), 1e14, rtol=1e-12)


@requires_file(athenapk_disk)
def test_load_cylindrical():
    # Load a cylindrical dataset of a full disk
    ds = data_dir_load(athenapk_disk)

    # Check that the domain edges match r in [0.5,2.0], theta in [0, 2pi]
    assert_equal(ds.domain_left_edge.in_units("code_length").v[:2], (0.5, 0))
    assert_equal(ds.domain_right_edge.in_units("code_length").v[:2], (2.0, 2 * np.pi))

"""
CF-Radial frontend tests



"""
import os
import shutil
import tempfile

from yt.frontends.cf_radial.data_structures import CFRadialDataset
from yt.testing import (
    assert_almost_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds

cf = "CfRadialGrid/grid1.nc"


@requires_file(cf)
def test_units_override():
    units_override_check(cf)


@requires_file(cf)
def test_CFRadialDataset():
    assert isinstance(data_dir_load(cf), CFRadialDataset)


_fields_xarray = ["reflectivity", "velocity", "gate_id", "differential_phase", "ROI"]


@requires_file(cf)
def test_xarray_fields():
    ds = data_dir_load(cf)
    for field in _fields_xarray:
        check_xarray_fields(ds, field)


def check_xarray_fields(ds, field):
    assert field in ds._handle.variables.keys()


_fields_cfradial = ["reflectivity", "velocity", "gate_id", "differential_phase", "ROI"]


@requires_file(cf)
def test_cfradial_fields():
    ds = data_dir_load(cf)
    for field in _fields_cfradial:
        check_fields(ds, field)


def check_fields(ds, field):
    assert ("cf_radial", field) in ds.field_list


_fields_units = {
    "reflectivity": "dimensionless",
    "velocity": "m/s",
    "differential_phase": "degree",
    "gate_id": "dimensionless",
    "ROI": "m",
}


@requires_file(cf)
def test_units():
    ds = data_dir_load(cf)
    ad = ds.all_data()
    for field in _fields_units.keys():
        check_field_units(ad, field, _fields_units[field])


def check_field_units(ad, field, value):
    assert str(ad["cf_radial", field].units) == value


@requires_file(cf)
def test_origin_latitude():
    ds = data_dir_load(cf)
    assert_almost_equal(ds.origin_latitude.values, 36.49120001)


@requires_file(cf)
def test_origin_longitude():
    ds = data_dir_load(cf)
    assert_almost_equal(ds.origin_longitude.values, -97.5939)


@requires_file(cf)
def test_domain_dimensions():
    ds = data_dir_load(cf)
    domain_dim_array = [251, 251, 46]
    assert_equal(ds.domain_dimensions, domain_dim_array)


@requires_file(cf)
def test_domain_center():
    ds = data_dir_load(cf)
    domain_center_array = [0.0, 0.0, 7500.0]
    assert_equal(ds.domain_center, domain_center_array)


@requires_file(cf)
def test_domain_left_edge():
    ds = data_dir_load(cf)
    domain_left_array = [-50000.0, -50000.0, 0.0]
    assert_equal(ds.domain_left_edge, domain_left_array)


@requires_file(cf)
def test_domain_right_edge():
    ds = data_dir_load(cf)
    domain_right_array = [50000.0, 50000.0, 15000.0]
    assert_equal(ds.domain_right_edge, domain_right_array)


cf_nongridded = "CfRadialGrid/swx_20120520_0641.nc"


@requires_ds(cf_nongridded)
def test_gridding():
    # loads up a radial dataset, which triggers the gridding.

    # create temporary directory and grid file
    tempdir = tempfile.mkdtemp()
    grid_file = os.path.join(tempdir, "temp_grid.nc")

    # this load will trigger the re-gridding and write out the gridded file
    # from which data will be loaded.
    ds = data_dir_load(cf_nongridded, kwargs={"storage_filename": grid_file})
    assert os.path.exists(grid_file)

    # check that the cartesian fields exist now
    for field in ["x", "y", "z"]:
        check_xarray_fields(ds, field)

    shutil.rmtree(tempdir)

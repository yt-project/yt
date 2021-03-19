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
    requires_module,
    units_override_check,
)
from yt.utilities.answer_testing.framework import data_dir_load, requires_ds

cf = "CfRadialGrid/grid1.nc"  # an already gridded cfradial file
cf_nongridded = (
    "CfRadialGrid/swx_20120520_0641.nc"  # cfradial file without cartesian grid
)
_fields_cfradial = ["reflectivity", "velocity", "gate_id", "differential_phase", "ROI"]

_fields_units = {
    "reflectivity": "dimensionless",
    "velocity": "m/s",
    "differential_phase": "degree",
    "gate_id": "dimensionless",
    "ROI": "m",
}


@requires_file(cf)
@requires_module("xarray")
def test_units_override():
    units_override_check(cf)


@requires_file(cf)
@requires_module("xarray")
def test_cf_radial_gridded():
    ds = data_dir_load(cf)
    assert isinstance(ds, CFRadialDataset)

    check_domain(ds)
    check_origin_latitude_longitude(ds)
    ad = ds.all_data()

    for field in _fields_cfradial:
        check_fields(ds, field)
        check_field_units(ad, field, _fields_units[field])


def check_fields(ds, field):
    assert ("cf_radial", field) in ds.field_list
    assert field in ds._handle.variables.keys()


def check_field_units(ad, field, value):
    assert str(ad["cf_radial", field].units) == value


def check_origin_latitude_longitude(ds):
    assert_almost_equal(ds.origin_latitude.values, 36.49120001)
    assert_almost_equal(ds.origin_longitude.values, -97.5939)


def check_domain(ds):
    domain_dim_array = [251, 251, 46]
    assert_equal(ds.domain_dimensions, domain_dim_array)

    domain_center_array = [0.0, 0.0, 7500.0]
    assert_equal(ds.domain_center, domain_center_array)

    domain_left_array = [-50000.0, -50000.0, 0.0]
    assert_equal(ds.domain_left_edge, domain_left_array)

    domain_right_array = [50000.0, 50000.0, 15000.0]
    assert_equal(ds.domain_right_edge, domain_right_array)


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
        assert field in ds._handle.variables.keys()

    shutil.rmtree(tempdir)

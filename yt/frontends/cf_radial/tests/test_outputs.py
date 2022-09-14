"""
CF-Radial frontend tests



"""
import os
import shutil
import tempfile

import numpy as np

from yt.frontends.cf_radial.data_structures import CFRadialDataset
from yt.testing import (
    assert_almost_equal,
    assert_equal,
    requires_file,
    requires_module,
    units_override_check,
)
from yt.utilities.answer_testing.framework import (
    data_dir_load,
    requires_ds,
    small_patch_amr,
)

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


@requires_module("xarray")
@requires_file(cf)
def test_units_override():
    units_override_check(cf)


@requires_module("xarray")
@requires_file(cf)
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
    with ds._handle() as xr_ds_handle:
        assert field in xr_ds_handle.variables.keys()


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


@requires_module("xarray")
@requires_file(cf_nongridded)
def test_auto_gridding():
    # loads up a radial dataset, which triggers the gridding.

    # create temporary directory and grid file
    tempdir = tempfile.mkdtemp()
    grid_file = os.path.join(tempdir, "temp_grid.nc")

    # this load will trigger the re-gridding and write out the gridded file
    # from which data will be loaded. With default gridding params, this takes
    # on the order of 10s, but since we are not testing actual output here, we
    # can decrease the resolution to speed it up.
    grid_shape = (10, 10, 10)
    ds = data_dir_load(
        cf_nongridded, kwargs={"storage_filename": grid_file, "grid_shape": grid_shape}
    )
    assert os.path.exists(grid_file)

    # check that the cartesian fields exist now
    with ds._handle() as xr_ds_handle:
        on_disk_fields = xr_ds_handle.variables.keys()
        for field in ["x", "y", "z"]:
            assert field in on_disk_fields

    assert all(ds.domain_dimensions == grid_shape)

    # check that we can load the gridded file too
    ds = data_dir_load(grid_file)
    assert isinstance(ds, CFRadialDataset)

    shutil.rmtree(tempdir)


@requires_module("xarray")
@requires_file(cf_nongridded)
def test_grid_parameters():
    # checks that the gridding parameters are used and that conflicts in parameters
    # are resolved as expected.
    tempdir = tempfile.mkdtemp()
    grid_file = os.path.join(tempdir, "temp_grid_params.nc")

    # check that the grid parameters work
    cfkwargs = {
        "storage_filename": grid_file,
        "grid_shape": (10, 10, 10),
        "grid_limit_x": (-10000, 10000),
        "grid_limit_y": (-10000, 10000),
        "grid_limit_z": (500, 20000),
    }
    ds = data_dir_load(cf_nongridded, kwargs=cfkwargs)

    expected_width = []
    for dim in "xyz":
        minval, maxval = cfkwargs[f"grid_limit_{dim}"]
        expected_width.append(maxval - minval)
    expected_width = np.array(expected_width)

    actual_width = ds.domain_width.to("m").value
    assert all(expected_width == actual_width)
    assert all(ds.domain_dimensions == cfkwargs["grid_shape"])

    # check the grid parameter conflicts

    # on re-load with default grid params it will reload storage_filename if
    # it exists. Just checking that this runs...
    _ = data_dir_load(cf_nongridded, kwargs={"storage_filename": grid_file})

    # if storage_filename exists, grid parameters are ignored (with a warning)
    # and the domain_dimensions will match the original
    new_kwargs = {"storage_filename": grid_file, "grid_shape": (15, 15, 15)}
    ds = data_dir_load(cf_nongridded, kwargs=new_kwargs)
    assert all(ds.domain_dimensions == cfkwargs["grid_shape"])

    # if we overwrite, the regridding should run and the dimensions should match
    # the desired dimensions
    new_kwargs["storage_overwrite"] = True
    ds = data_dir_load(cf_nongridded, kwargs=new_kwargs)
    assert all(ds.domain_dimensions == new_kwargs["grid_shape"])

    shutil.rmtree(tempdir)


@requires_module("xarray")
@requires_ds(cf)
def test_cfradial_grid_field_values():
    ds = data_dir_load(cf)
    fields_to_check = [("cf_radial", field) for field in _fields_cfradial]
    wtfield = ("cf_radial", "reflectivity")
    for test in small_patch_amr(
        ds, fields_to_check, input_center=ds.domain_center, input_weight=wtfield
    ):
        test_cfradial_grid_field_values.__name__ = test.description
        yield test

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
from yt.utilities.answer_testing.framework import data_dir_load

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
    for field in ["x", "y", "z"]:
        assert field in ds._handle.variables.keys()

    assert all(ds.domain_dimensions == grid_shape)

    # check that we can load the gridded file too
    ds = data_dir_load(grid_file)
    assert isinstance(ds, CFRadialDataset)

    shutil.rmtree(tempdir)


@requires_file(cf_nongridded)
def test_grid_parameters():
    # check that the ds.domain matches what we give for grid limits
    tempdir = tempfile.mkdtemp()
    grid_file = os.path.join(tempdir, "temp_grid_params.nc")

    # check that the other grid parameters work
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

    shutil.rmtree(tempdir)


@requires_file(cf_nongridded)
def test_grid_param_conflicts():

    # create temporary directory and grid file
    tempdir = tempfile.mkdtemp()
    grid_file = os.path.join(tempdir, "temp_grid_conflicts.nc")
    grid_shape = (10, 10, 10)
    ds = data_dir_load(
        cf_nongridded, kwargs={"storage_filename": grid_file, "grid_shape": grid_shape}
    )
    assert os.path.exists(grid_file)

    # on re-load with default grid params it will reload storage_filename if
    # it exists.
    _ = data_dir_load(cf_nongridded, kwargs={"storage_filename": grid_file})

    # if storage_filename exists, grid parameters are ignored (with a warning)
    # and the domain_dimensions will match the original
    bad_kwargs = {"storage_filename": grid_file, "grid_shape": (15, 15, 15)}
    ds = data_dir_load(cf_nongridded, kwargs=bad_kwargs)
    assert all(ds.domain_dimensions == grid_shape)
    ds._handle.close()  # TODO: handle is staying open. need to fix that.

    # if we overwrite, the regridding should run and the dimensions should match
    # the desired dimensions

    cfkwargs = {
        "storage_filename": grid_file,
        "grid_shape": (15, 15, 15),
        "storage_overwrite": True,
    }
    ds = data_dir_load(cf_nongridded, kwargs=cfkwargs)
    print(ds.domain_dimensions)
    print(cfkwargs["grid_shape"])
    assert all(ds.domain_dimensions == cfkwargs["grid_shape"])

    shutil.rmtree(tempdir)

"""
CF-Radial frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    assert_almost_equal, \
    assert_equal, \
    requires_file, \
    units_override_check
from yt.utilities.answer_testing.framework import \
    small_patch_amr, \
    data_dir_load
from ..data_structures import CFRadialDataset


cf = "CfRadialGrid/grid1.nc"

@requires_file(cf)
def test_units_override():
    units_override_check(cf)

@requires_file(cf)
def test_CFRadialDataset():
    assert isinstance(data_dir_load(cf), CFRadialDataset)

_fields_xarray = ["reflectivity", "velocity", "gate_id",
                  "differential_phase", "ROI"]
@requires_file(cf)
def test_xarray_fields():
    ds = data_dir_load(cf)
    for field in _fields_xarray:
        check_xarray_fields(ds, field)

def check_xarray_fields(ds, field):
    assert field in ds._handle.variables.keys()

_fields_cfradial = ["reflectivity", "velocity", "gate_id",
                    "differential_phase", "ROI"]
@requires_file(cf)
def test_cfradial_fields():
    ds = data_dir_load(cf)
    for field in _fields_cfradial:
        check_fields(ds, field)

def check_fields(ds, field):
    assert ("cf_radial", field) in ds.field_list

_fields_units = {"reflectivity": "dBZ", "velocity": "m/s",
                 "differential_phase": "degree",
                 "gate_id": "dimensionless", "ROI": "dimensionless"}
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
def test_domain_dimensions()::

@requires_files(cf)
def test_domain_left_edge():

@requires_files(cf)
def test_domain_right_edge():

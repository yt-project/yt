"""
Swift frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt import load

from yt.testing import requires_file
from yt.frontends.swift.api import SwiftDataset
from yt.testing import assert_almost_equal
from yt.utilities.on_demand_imports import _h5py as h5py


keplerian_ring = "KeplerianRing/keplerian_ring_test_data.hdf5"
EAGLE_6 = "EAGLE_6/eagle_00010.hdf5"

# Combined the tests for loading a file and ensuring the units have been
# implemented correctly to save time on re-loading a dataset
@requires_file(keplerian_ring)
def test_non_cosmo_dataset():
    ds = load(keplerian_ring)
    assert(type(ds) == SwiftDataset)

    # now we have the file, lets test the units
    field = ('gas', 'density')
    ad = ds.all_data()
    ytval = ad[field]

    # load some data the old fashioned way
    fh = h5py.File(keplerian_ring, "r")
    part_data = fh['PartType0']
    units = fh["Units"]
    units = dict(units.attrs)
    conv_factor = float(units["Unit mass in cgs (U_M)"])
    conv_factor /= float(units["Unit length in cgs (U_L)"])**3
    h5val = part_data["Density"][:].astype("float64") * conv_factor
    fh.close()

    # make sure we are comparing fair units
    assert(str(ytval.units) == 'g/cm**3')

    # make sure the actual values are the same
    assert_almost_equal(h5val, ytval.d)

    return

@requires_file(EAGLE_6)
def test_cosmo_dataset():
    ds = load(EAGLE_6)
    assert(type(ds) == SwiftDataset)

    # lets test the cosmology units
    field = ('gas', 'density')
    ad = ds.all_data()
    ytval = ad[field]

    # load some data the old fashioned way
    fh = h5py.File(keplerian_ring, "r")
    part_data = fh['PartType0']
    units = fh["Units"]
    units = dict(units.attrs)
    conv_factor = float(units["Unit mass in cgs (U_M)"])
    conv_factor /= float(units["Unit length in cgs (U_L)"])**3
    h5val = part_data["Density"][:].astype("float64") * conv_factor
    fh.close()

    # make sure we are comparing fair units
    assert(str(ytval.units) == 'g/cm**3')

    # make sure the actual values are the same
    assert_almost_equal(h5val, ytval.d)

    return

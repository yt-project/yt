"""
GadgetFOF frontend tests using gadget_fof datasets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.frontends.gadget_fof.api import \
    GadgetFOFDataset
from yt.testing import \
    requires_file, \
    assert_equal, \
    assert_array_equal
from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load

_fields = ("particle_position_x", "particle_position_y",
           "particle_position_z", "particle_velocity_x",
           "particle_velocity_y", "particle_velocity_z",
           "particle_mass", "particle_identifier")

# a dataset with empty files
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"


@requires_ds(g5)
def test_fields_g5():
    for field in _fields:
        yield FieldValuesTest(g5, field, particle_type=True)


@requires_ds(g42)
def test_fields_g42():
    for field in _fields:
        yield FieldValuesTest(g42, field, particle_type=True)

@requires_file(g42)
def test_GadgetFOFDataset():
    assert isinstance(data_dir_load(g42), GadgetFOFDataset)

# fof/subhalo catalog with member particles
g298 = "gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5"

@requires_file(g298)
def test_subhalos():
    ds = data_dir_load(g298)
    total_sub = 0
    total_int = 0
    for hid in range(0, ds.index.particle_count["Group"]):
        my_h = ds.halo("Group", hid)
        h_ids = my_h["ID"]
        for sid in range(my_h["subhalo_number"]):
            my_s = ds.halo("Subhalo", (my_h.particle_identifier, sid))
            total_sub += my_s["ID"].size
            total_int += np.intersect1d(h_ids, my_s["ID"]).size

    # Test that all subhalo particles are contained within
    # their parent group.
    yield assert_equal, total_sub, total_int

@requires_file(g298)
def test_halomasses():
    ds = data_dir_load(g298)
    ad = ds.all_data()
    for ptype in ["Group", "Subhalo"]:
        nhalos = ds.index.particle_count[ptype]
        mass = ds.arr(np.zeros(nhalos), "code_mass")
        for i in range(nhalos):
            halo = ds.halo(ptype, i)
            mass[i] = halo.mass

        # Check that masses from halo containers are the same
        # as the array of all masses.  This will test getting
        # scalar fields for halos correctly.
        yield assert_array_equal, ad[ptype, "particle_mass"], mass

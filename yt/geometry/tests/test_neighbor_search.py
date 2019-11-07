"""
Tests for neighbor finding



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.fields.particle_fields import \
    add_nearest_neighbor_field, \
    add_nearest_neighbor_value_field
from yt.testing import \
    fake_particle_ds, \
    assert_equal, \
    assert_array_almost_equal


def test_neighbor_search():
    np.random.seed(0x4d3d3d3)
    ds = fake_particle_ds(npart = 16**3)
    ds.periodicity = (True, True, True)
    ds.index
    fn, = add_nearest_neighbor_field("all", "particle_position", ds)
    dd = ds.all_data()
    nearest_neighbors = dd[fn]
    pos = dd["particle_position"]
    all_neighbors = np.zeros_like(nearest_neighbors)
    any_eq = np.zeros(pos.shape[0], dtype='bool')
    min_in = np.zeros(pos.shape[0], dtype='int64')
    for i in range(pos.shape[0]):
        dd.set_field_parameter("center", pos[i,:])
        #radius = dd["particle_radius"]
        #radius.sort()
        r2 = (pos[:,0]*pos[:,0])*0
        for j in range(3):
            DR = (pos[i,j] - pos[:,j])
            DRo = DR.copy()
            DR[DRo >  ds.domain_width[j]/2.0] -= ds.domain_width[j]
            DR[DRo < -ds.domain_width[j]/2.0] += ds.domain_width[j]
            r2 += DR*DR
        radius = np.sqrt(r2)
        radius.sort()
        assert(radius[0] == 0.0)
        all_neighbors[i] = radius[63]
        any_eq[i] = np.any( np.abs(radius - nearest_neighbors[i]) < 1e-7 )
        min_in[i] = np.argmin(np.abs(radius - nearest_neighbors[i]))
        #if i == 34: raise RuntimeError
        #dd.field_data.pop(("all", "particle_radius"))
    assert_equal((min_in == 63).sum(), min_in.size)
    assert_array_almost_equal(nearest_neighbors, all_neighbors)

def test_neighbor_value_search():
    np.random.seed(0x4d3d3d3)
    ds = fake_particle_ds(npart = 16**3, over_refine_factor = 2)
    ds.periodicity = (True, True, True)
    ds.index
    fn, = add_nearest_neighbor_value_field("all", "particle_position", "particle_mass",
                                           ds.field_info)
    dd = ds.all_data()
    # Set up our positions onto which the field will be deposited
    index_pos = np.array([dd["index",ax] for ax in 'xyz']) * dd["index","x"].uq
    particle_pos = dd["particle_position"]
    values_in = dd["particle_mass"]
    values_out = dd[fn]
    for i in range(index_pos.shape[0]):
        r2 = particle_pos[:,0]*0
        r2 = r2 * r2
        for j in range(3):
            DR = (index_pos[i,j] - particle_pos[:,j])
            DRo = DR.copy()
            DR[DRo >  ds.domain_width[j]/2.0] -= ds.domain_width[j]
            DR[DRo < -ds.domain_width[j]/2.0] += ds.domain_width[j]
            r2 += DR*DR
        radius = np.sqrt(r2)
        assert(values_in[np.argmin(radius)] == values_out[i])

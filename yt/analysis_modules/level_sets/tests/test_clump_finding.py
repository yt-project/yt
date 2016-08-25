"""
Clump finder tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.analysis_modules.level_sets.api import \
    Clump, \
    find_clumps, \
    get_lowest_clumps
from yt.frontends.stream.api import \
    load_uniform_grid
from yt.testing import \
    assert_array_equal, \
    assert_equal

def test_clump_finding():
    n_c = 8
    n_p = 1
    dims = (n_c, n_c, n_c)

    density = np.ones(dims)
    high_rho = 10.
    # add a couple disconnected density enhancements
    density[2, 2, 2] = high_rho
    density[6, 6, 6] = high_rho

    # put a particle at the center of one of them
    dx = 1. / n_c
    px = 2.5 * dx * np.ones(n_p)
    
    data = {"density": density,
            "particle_mass": np.ones(n_p),
            "particle_position_x": px,
            "particle_position_y": px,
            "particle_position_z": px,
            "number_of_particles": n_p}

    ds = load_uniform_grid(data, dims)

    ad = ds.all_data()
    master_clump = Clump(ad, ("gas", "density"))
    master_clump.add_validator("min_cells", 1)

    find_clumps(master_clump, 0.5, 2. * high_rho, 10.)

    # there should be two children
    assert_equal(len(master_clump.children), 2)

    leaf_clumps = get_lowest_clumps(master_clump)
    # two leaf clumps
    assert_equal(len(leaf_clumps), 2)


    # check some clump fields
    assert_equal(master_clump.children[0]["density"][0].size, 1)
    assert_equal(master_clump.children[0]["density"][0], ad["density"].max())
    assert_equal(master_clump.children[0]["particle_mass"].size, 1)
    assert_array_equal(master_clump.children[0]["particle_mass"], ad["particle_mass"])
    assert_equal(master_clump.children[1]["density"][0].size, 1)
    assert_equal(master_clump.children[1]["density"][0], ad["density"].max())
    assert_equal(master_clump.children[1]["particle_mass"].size, 0)

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
import os
import shutil
import tempfile


from yt.analysis_modules.level_sets.api import \
    Clump, \
    find_clumps, \
    get_lowest_clumps
from yt.convenience import \
    load
from yt.frontends.stream.api import \
    load_uniform_grid
from yt.testing import \
    assert_array_equal, \
    assert_equal, \
    requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load

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

i30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
@requires_file(i30)
def test_clump_tree_save():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = data_dir_load(i30)
    data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.],
                          (8, 'kpc'), (1, 'kpc'))

    field = ("gas", "density")
    step = 2.0
    c_min = 10**np.floor(np.log10(data_source[field]).min()  )
    c_max = 10**np.floor(np.log10(data_source[field]).max()+1)

    master_clump = Clump(data_source, field)
    master_clump.add_info_item("center_of_mass")
    master_clump.add_validator("min_cells", 20)

    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = get_lowest_clumps(master_clump)

    fn = master_clump.save_as_dataset(fields=["density", "x", "y", "z",
                                              "particle_mass"])
    ds2 = load(fn)

    # compare clumps in the tree
    t1 = [c for c in master_clump]
    t2 = [c for c in ds2.tree]
    mt1 = ds.arr([c.info["cell_mass"][1] for c in t1])
    mt2 = ds2.arr([c["clump", "cell_mass"] for c in t2])
    it1 = np.argsort(mt1).d.astype(int)
    it2 = np.argsort(mt2).d.astype(int)
    assert_array_equal(mt1[it1], mt2[it2])

    for i1, i2 in zip(it1, it2):
        ct1 = t1[i1]
        ct2 = t2[i2]
        assert_array_equal(ct1["gas", "density"],
                           ct2["grid", "density"])
        assert_array_equal(ct1["all", "particle_mass"],
                           ct2["all", "particle_mass"])

    # compare leaf clumps
    c1 = [c for c in leaf_clumps]
    c2 = [c for c in ds2.leaves]
    mc1 = ds.arr([c.info["cell_mass"][1] for c in c1])
    mc2 = ds2.arr([c["clump", "cell_mass"] for c in c2])
    ic1 = np.argsort(mc1).d.astype(int)
    ic2 = np.argsort(mc2).d.astype(int)
    assert_array_equal(mc1[ic1], mc2[ic2])

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

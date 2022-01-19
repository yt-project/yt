import os
import shutil
import tempfile

import numpy as np

from yt.data_objects.level_sets.api import Clump, add_clump_info, find_clumps
from yt.data_objects.level_sets.clump_info_items import clump_info_registry
from yt.fields.derived_field import ValidateParameter
from yt.loaders import load, load_uniform_grid
from yt.testing import assert_array_equal, assert_equal, requires_file
from yt.utilities.answer_testing.framework import data_dir_load


def test_clump_finding():
    n_c = 8
    n_p = 1
    dims = (n_c, n_c, n_c)

    density = np.ones(dims)
    high_rho = 10.0
    # add a couple disconnected density enhancements
    density[2, 2, 2] = high_rho
    density[6, 6, 6] = high_rho

    # put a particle at the center of one of them
    dx = 1.0 / n_c
    px = 2.5 * dx * np.ones(n_p)

    data = {
        "density": density,
        "particle_mass": np.ones(n_p),
        "particle_position_x": px,
        "particle_position_y": px,
        "particle_position_z": px,
    }

    ds = load_uniform_grid(data, dims)

    ad = ds.all_data()
    master_clump = Clump(ad, ("gas", "density"))
    master_clump.add_validator("min_cells", 1)

    def _total_volume(clump):
        total_vol = clump.data.quantities.total_quantity(
            [("index", "cell_volume")]
        ).in_units("cm**3")
        return "Cell Volume: %6e cm**3.", total_vol

    add_clump_info("total_volume", _total_volume)
    master_clump.add_info_item("total_volume")

    find_clumps(master_clump, 0.5, 2.0 * high_rho, 10.0)

    # there should be two children
    assert_equal(len(master_clump.children), 2)

    leaf_clumps = master_clump.leaves

    for l in leaf_clumps:
        keys = l.info.keys()
        assert "total_cells" in keys
        assert "cell_mass" in keys
        assert "max_grid_level" in keys
        assert "total_volume" in keys

    # two leaf clumps
    assert_equal(len(leaf_clumps), 2)

    # check some clump fields
    assert_equal(master_clump.children[0][("gas", "density")][0].size, 1)
    assert_equal(
        master_clump.children[0][("gas", "density")][0], ad[("gas", "density")].max()
    )
    assert_equal(master_clump.children[0][("all", "particle_mass")].size, 1)
    assert_array_equal(
        master_clump.children[0][("all", "particle_mass")], ad[("all", "particle_mass")]
    )
    assert_equal(master_clump.children[1][("gas", "density")][0].size, 1)
    assert_equal(
        master_clump.children[1][("gas", "density")][0], ad[("gas", "density")].max()
    )
    assert_equal(master_clump.children[1][("all", "particle_mass")].size, 0)

    # clean up global registry to avoid polluting other tests
    del clump_info_registry["total_volume"]


i30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(i30)
def test_clump_tree_save():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    ds = data_dir_load(i30)
    data_source = ds.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc"))

    field = ("gas", "density")
    step = 2.0
    c_min = 10 ** np.floor(np.log10(data_source[field]).min())
    c_max = 10 ** np.floor(np.log10(data_source[field]).max() + 1)

    master_clump = Clump(data_source, field)
    master_clump.add_info_item("center_of_mass")
    master_clump.add_validator("min_cells", 20)

    find_clumps(master_clump, c_min, c_max, step)
    leaf_clumps = master_clump.leaves

    fn = master_clump.save_as_dataset(
        fields=[
            ("gas", "density"),
            ("index", "x"),
            ("index", "y"),
            ("index", "z"),
            ("all", "particle_mass"),
        ]
    )
    ds2 = load(fn)

    # compare clumps in the tree
    t1 = [c for c in master_clump]
    t2 = [c for c in ds2.tree]
    mt1 = ds.arr([c.info["cell_mass"][1] for c in t1])
    mt2 = ds2.arr([c["clump", "cell_mass"] for c in t2])
    it1 = np.array(np.argsort(mt1).astype(int))
    it2 = np.array(np.argsort(mt2).astype(int))
    assert_array_equal(mt1[it1], mt2[it2])

    for i1, i2 in zip(it1, it2):
        ct1 = t1[i1]
        ct2 = t2[i2]
        assert_array_equal(ct1["gas", "density"], ct2["grid", "density"])
        assert_array_equal(ct1["all", "particle_mass"], ct2["all", "particle_mass"])

    # compare leaf clumps
    c1 = [c for c in leaf_clumps]
    c2 = [c for c in ds2.leaves]
    mc1 = ds.arr([c.info["cell_mass"][1] for c in c1])
    mc2 = ds2.arr([c["clump", "cell_mass"] for c in c2])
    ic1 = np.array(np.argsort(mc1).astype(int))
    ic2 = np.array(np.argsort(mc2).astype(int))
    assert_array_equal(mc1[ic1], mc2[ic2])

    os.chdir(curdir)
    shutil.rmtree(tmpdir)


i30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(i30)
def test_clump_field_parameters():
    """
    Make sure clump finding on fields with field parameters works.
    """

    def _also_density(field, data):
        factor = data.get_field_parameter("factor")
        return factor * data[("gas", "density")]

    ds = data_dir_load(i30)
    ds.add_field(
        ("gas", "also_density"),
        function=_also_density,
        units=ds.fields.gas.density.units,
        sampling_type="cell",
        validators=[ValidateParameter("factor")],
    )
    data_source = ds.disk([0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc"))
    data_source.set_field_parameter("factor", 1)

    step = 2.0
    field = ("gas", "density")
    c_min = 10 ** np.floor(np.log10(data_source[field]).min())
    c_max = 10 ** np.floor(np.log10(data_source[field]).max() + 1)

    master_clump_1 = Clump(data_source, ("gas", "density"))
    master_clump_1.add_validator("min_cells", 20)
    master_clump_2 = Clump(data_source, ("gas", "also_density"))
    master_clump_2.add_validator("min_cells", 20)

    find_clumps(master_clump_1, c_min, c_max, step)
    find_clumps(master_clump_2, c_min, c_max, step)
    leaf_clumps_1 = master_clump_1.leaves
    leaf_clumps_2 = master_clump_2.leaves

    for c1, c2 in zip(leaf_clumps_1, leaf_clumps_2):
        assert_array_equal(c1["gas", "density"], c2["gas", "density"])

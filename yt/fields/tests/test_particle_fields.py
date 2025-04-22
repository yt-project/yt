import numpy as np

from yt.testing import assert_allclose_units, requires_file, requires_module
from yt.utilities.answer_testing.framework import data_dir_load

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_module("h5py")
@requires_file(g30)
def test_relative_particle_fields():
    ds = data_dir_load(g30)
    offset = ds.arr([0.1, -0.2, 0.3], "code_length")
    c = ds.domain_center + offset
    sp = ds.sphere(c, (10, "kpc"))
    bv = ds.arr([1.0, 2.0, 3.0], "code_velocity")
    sp.set_field_parameter("bulk_velocity", bv)
    assert_allclose_units(
        sp["all", "relative_particle_position"], sp["all", "particle_position"] - c
    )
    assert_allclose_units(
        sp["all", "relative_particle_velocity"], sp["all", "particle_velocity"] - bv
    )


@requires_module("h5py")
@requires_file(g30)
def test_los_particle_fields():
    ds = data_dir_load(g30)
    offset = ds.arr([0.1, -0.2, 0.3], "code_length")
    c = ds.domain_center + offset
    sp = ds.sphere(c, (10, "kpc"))
    bv = ds.arr([1.0, 2.0, 3.0], "code_velocity")
    sp.set_field_parameter("bulk_velocity", bv)
    ax = [0.1, 0.2, -0.3]
    sp.set_field_parameter("axis", ax)
    ax /= np.linalg.norm(ax)
    vlos = (
        sp["all", "relative_particle_velocity_x"] * ax[0]
        + sp["all", "relative_particle_velocity_y"] * ax[1]
        + sp["all", "relative_particle_velocity_z"] * ax[2]
    )
    assert_allclose_units(sp["all", "particle_velocity_los"], vlos)
    sp.clear_data()
    ax = 2
    sp.set_field_parameter("axis", ax)
    vlos = sp["all", "relative_particle_velocity_z"]
    assert_allclose_units(sp["all", "particle_velocity_los"], vlos)

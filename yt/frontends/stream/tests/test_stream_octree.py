import numpy as np
import pytest
from numpy.testing import assert_equal

import yt

OCT_MASK_LIST = [
    8,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    8,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]


def test_octree():
    # See Issue #1272
    octree_mask = np.array(OCT_MASK_LIST, dtype=np.uint8)

    quantities = {}
    rng = np.random.default_rng(12345)
    quantities["gas", "density"] = rng.random((22, 1))
    quantities["gas", "dinos"] = (np.linspace(0, 1, 22).reshape(22, 1), "Msun")
    quantities["gas", "turtle"] = (lambda: np.linspace(0, 1, 22).reshape(22, 1), "Msun")

    bbox = np.array([[-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0]])

    ds = yt.load_octree(
        octree_mask=octree_mask,
        data=quantities,
        bbox=bbox,
        num_zones=1,
        partial_coverage=0,
    )

    proj = ds.proj(("gas", "density"), "x")
    proj["gas", "density"]

    assert_equal(ds.r[:]["ones"].size, 22)
    rho1 = quantities["gas", "density"].ravel()
    rho2 = ds.r[:]["density"].copy()
    rho1.sort()
    rho2.sort()
    assert_equal(rho1, rho2)

    # Make sure units aren't stripped
    assert str(ds.r[:]["dinos"].units) == "Msun"
    assert ds.r[:]["dinos"].min().to("Msun") == 0
    assert ds.r[:]["dinos"].max().to("Msun") == 1

    assert str(ds.r[:]["turtle"].units) == "Msun"
    assert ds.r[:]["turtle"].min().to("Msun") == 0
    assert ds.r[:]["turtle"].max().to("Msun") == 1


@pytest.mark.parametrize("build_pos", ["component_by_component", "combined"])
def test_octree_particles(build_pos):
    octree_mask = np.array(OCT_MASK_LIST, dtype=np.uint8)

    L = [-10, -10, -10]
    R = [10, 10, 10]

    quantities = {}
    rng = np.random.default_rng(12345)

    pos = rng.uniform(-10, 10, size=(100, 3))
    quantities["gas", "density"] = rng.random((22, 1))
    if build_pos == "component_by_component":
        quantities["io", "particle_position_x"] = pos[:, 0]
        quantities["io", "particle_position_y"] = pos[:, 1]
        quantities["io", "particle_position_z"] = pos[:, 2]
    else:
        quantities["io", "particle_position"] = pos
    quantities["io", "particle_mass"] = (np.linspace(0, 1, 100), "Msun")

    bbox = np.array([L, R]).T

    ds = yt.load_octree(
        octree_mask=octree_mask,
        data=quantities,
        bbox=bbox,
        num_zones=1,
        partial_coverage=0,
    )

    # Ensure "io" is registered as a particle type
    assert "io" in ds.particle_types
    assert "io" in ds.particle_types_raw

    # Ensure we can access the particle positions
    _pos = ds.r[:]["particle_position"]
    _pos_x = ds.r[:]["particle_position_x"]
    _pos_y = ds.r[:]["particle_position_y"]
    _pos_z = ds.r[:]["particle_position_z"]

    # Make sure other fields are defined
    mass = ds.r[:]["particle_mass"].to("Msun")
    assert mass.min() == 0
    assert mass.max() == 1

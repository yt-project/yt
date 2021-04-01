import numpy as np
from nose.tools import assert_raises

from yt.testing import (
    assert_almost_equal,
    assert_equal,
    fake_amr_ds,
    fake_particle_ds,
    fake_random_ds,
    requires_file,
)
from yt.utilities.answer_testing.framework import data_dir_load
from yt.utilities.exceptions import YTDimensionalityError
from yt.visualization.line_plot import LineBuffer

# This will test the "dataset access" method.


def test_box_creation():
    ds = fake_random_ds(32, length_unit=2)
    left_edge = ds.arr([0.2, 0.2, 0.2], "cm")
    right_edge = ds.arr([0.6, 0.6, 0.6], "cm")
    center = (left_edge + right_edge) / 2

    boxes = [
        ds.box(left_edge, right_edge),
        ds.box(0.5 * np.array(left_edge), 0.5 * np.array(right_edge)),
        ds.box((0.5 * left_edge).tolist(), (0.5 * right_edge).tolist()),
    ]

    region = ds.region(center, left_edge, right_edge)

    for b in boxes:
        assert_almost_equal(b.left_edge, region.left_edge)
        assert_almost_equal(b.right_edge, region.right_edge)
        assert_almost_equal(b.center, region.center)


def test_region_from_d():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    # We'll do a couple here

    # First, no string units
    reg1 = ds.r[0.2:0.3, 0.4:0.6, :]
    reg2 = ds.region([0.25, 0.5, 0.5], [0.2, 0.4, 0.0], [0.3, 0.6, 1.0])
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])

    # Now, string units in some -- 1.0 == cm
    reg1 = ds.r[(0.1, "cm"):(0.5, "cm"), :, (0.25, "cm"):(0.35, "cm")]
    reg2 = ds.region([0.3, 0.5, 0.3], [0.1, 0.0, 0.25], [0.5, 1.0, 0.35])
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])

    # Now, string units in some -- 1.0 == cm
    reg1 = ds.r[(0.1, "cm"):(0.5, "cm"), :, 0.25:0.35]
    reg2 = ds.region([0.3, 0.5, 0.3], [0.1, 0.0, 0.25], [0.5, 1.0, 0.35])
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])

    # And, lots of : usage!
    reg1 = ds.r[:, :, :]
    reg2 = ds.all_data()
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])

    # Test slice as an index
    reg1 = ds.r[0.1:0.8]
    reg2 = ds.region([0.45, 0.45, 0.45], [0.1, 0.1, 0.1], [0.8, 0.8, 0.8])
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])

    # Test with bad boundary initialization
    with assert_raises(RuntimeError):
        ds.r[0.3:0.1, 0.4:0.6, :]

    # Test region by creating an arbitrary grid
    reg1 = ds.r[0.15:0.55:16j, 0.25:0.65:32j, 0.35:0.75:64j]
    left_edge = np.array([0.15, 0.25, 0.35])
    right_edge = np.array([0.55, 0.65, 0.75])
    dims = np.array([16.0, 32.0, 64.0])
    reg2 = ds.arbitrary_grid(left_edge, right_edge, dims)
    assert_equal(reg1[("gas", "density")], reg2[("gas", "density")])


def test_accessing_all_data():
    # This will test first that we can access all_data, and next that we can
    # access it multiple times and get the *same object*.
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    dd = ds.all_data()
    assert_equal(ds.r[("gas", "density")], dd[("gas", "density")])
    # Now let's assert that it's the same object
    rho = ds.r[("gas", "density")]
    rho *= 2.0
    assert_equal(dd[("gas", "density")] * 2.0, ds.r[("gas", "density")])
    assert_equal(dd["gas", "density"] * 2.0, ds.r["gas", "density"])


def test_slice_from_r():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    sl1 = ds.r[0.5, :, :]
    sl2 = ds.slice("x", 0.5)
    assert_equal(sl1[("gas", "density")], sl2[("gas", "density")])

    frb1 = sl1.to_frb(width=1.0, height=1.0, resolution=(1024, 512))
    frb2 = ds.r[0.5, ::1024j, ::512j]
    assert_equal(frb1[("gas", "density")], frb2[("gas", "density")])

    # Test slice which doesn't cover the whole domain
    box = ds.box([0.0, 0.25, 0.25], [1.0, 0.75, 0.75])

    sl3 = ds.r[0.5, 0.25:0.75, 0.25:0.75]
    sl4 = ds.slice("x", 0.5, data_source=box)
    assert_equal(sl3[("gas", "density")], sl4[("gas", "density")])

    frb3 = sl3.to_frb(width=0.5, height=0.5, resolution=(1024, 512))
    frb4 = ds.r[0.5, 0.25:0.75:1024j, 0.25:0.75:512j]
    assert_equal(frb3[("gas", "density")], frb4[("gas", "density")])


def test_point_from_r():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    pt1 = ds.r[0.5, 0.3, 0.1]
    pt2 = ds.point([0.5, 0.3, 0.1])
    assert_equal(pt1[("gas", "density")], pt2[("gas", "density")])

    # Test YTDimensionalityError
    with assert_raises(YTDimensionalityError) as ex:
        ds.r[0.5, 0.1]
    assert_equal(str(ex.exception), "Dimensionality specified was 2 but we need 3")


def test_ray_from_r():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    ray1 = ds.r[(0.1, 0.2, 0.3):(0.4, 0.5, 0.6)]
    ray2 = ds.ray((0.1, 0.2, 0.3), (0.4, 0.5, 0.6))
    assert_equal(ray1[("gas", "density")], ray2[("gas", "density")])

    ray3 = ds.r[0.5 * ds.domain_left_edge : 0.5 * ds.domain_right_edge]
    ray4 = ds.ray(0.5 * ds.domain_left_edge, 0.5 * ds.domain_right_edge)
    assert_equal(ray3[("gas", "density")], ray4[("gas", "density")])

    start = [(0.1, "cm"), 0.2, (0.3, "cm")]
    end = [(0.5, "cm"), (0.4, "cm"), 0.6]
    ray5 = ds.r[start:end]
    start_arr = [ds.quan(0.1, "cm"), ds.quan(0.2, "cm"), ds.quan(0.3, "cm")]
    end_arr = [ds.quan(0.5, "cm"), ds.quan(0.4, "cm"), ds.quan(0.6, "cm")]
    ray6 = ds.ray(start_arr, end_arr)
    assert_equal(ray5[("gas", "density")], ray6[("gas", "density")])

    ray7 = ds.r[start:end:500j]
    ray8 = LineBuffer(ds, [0.1, 0.2, 0.3], [0.5, 0.4, 0.6], 500)
    assert_equal(ray7[("gas", "density")], ray8[("gas", "density")])


def test_ortho_ray_from_r():
    ds = fake_amr_ds(fields=["density"], units=["g/cm**3"])
    ray1 = ds.r[:, 0.3, 0.2]
    ray2 = ds.ortho_ray("x", [0.3, 0.2])
    assert_equal(ray1[("gas", "density")], ray2[("gas", "density")])

    # the y-coord is funny so test it too
    ray3 = ds.r[0.3, :, 0.2]
    ray4 = ds.ortho_ray("y", [0.2, 0.3])
    assert_equal(ray3[("gas", "density")], ray4[("gas", "density")])

    # Test ray which doesn't cover the whole domain
    box = ds.box([0.25, 0.0, 0.0], [0.75, 1.0, 1.0])
    ray5 = ds.r[0.25:0.75, 0.3, 0.2]
    ray6 = ds.ortho_ray("x", [0.3, 0.2], data_source=box)
    assert_equal(ray5[("gas", "density")], ray6[("gas", "density")])

    # Test fixed-resolution rays
    ray7 = ds.r[0.25:0.75:100j, 0.3, 0.2]
    ray8 = LineBuffer(ds, [0.2525, 0.3, 0.2], [0.7475, 0.3, 0.2], 100)
    assert_equal(ray7[("gas", "density")], ray8[("gas", "density")])


def test_particle_counts():
    ds = fake_random_ds(16, particles=100)
    assert ds.particle_type_counts == {"io": 100}

    pds = fake_particle_ds(npart=128)
    assert pds.particle_type_counts == {"io": 128}


g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(g30)
def test_checksum():
    assert fake_random_ds(16).checksum == "notafile"
    assert data_dir_load(g30).checksum == "6169536e4b9f737ce3d3ad440df44c58"

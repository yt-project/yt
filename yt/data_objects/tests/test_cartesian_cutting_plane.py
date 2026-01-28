import itertools

import matplotlib.pyplot as plt
import numpy as np
import pytest
import unyt

import yt
from yt.geometry.coordinates.spherical_coordinates import spherical_to_cartesian
from yt.testing import fake_amr_ds


def test_cartesian_cutting_plane():
    ds = fake_amr_ds(geometry="spherical")
    normal = np.array([0.0, 0.0, 1.0])
    plane_center = np.array([0.0, 0.0, 0.5])
    slc = ds.cartesian_cutting(normal, plane_center)
    frb = slc.to_frb(2.0, 800)
    bvals = frb[("index", "r")]
    mask = frb.get_mask(("index", "r"))
    # note: the min value of r on the plane will be the z value of the
    # plane center. how close it is to the correct answer will depend
    # on the size of the elements.
    assert np.allclose(bvals[mask].min().d, plane_center[2], atol=0.02)


def _get_spherical_uniform_grid(shp, bbox, axis_order):

    data = {"density": np.random.random(shp)}

    def _z(field, data):
        r = data["index", "r"]
        theta = data["index", "theta"]
        phi = data["index", "phi"]
        _, _, z = spherical_to_cartesian(r, theta, phi)
        return unyt.unyt_array(z, r.units)

    ds = yt.load_uniform_grid(
        data,
        shp,
        bbox=bbox,
        geometry="spherical",
        axis_order=axis_order,
        length_unit="m",
    )

    ds.add_field(
        name=("index", "z_val"), function=_z, sampling_type="cell", take_log=False
    )
    return ds


@pytest.fixture
def spherical_ds():

    shp = (32, 32, 32)
    bbox = np.array([[0.0, 1.0], [0, np.pi], [0, 2 * np.pi]])
    ax_order = ("r", "theta", "phi")
    return _get_spherical_uniform_grid(shp, bbox, ax_order)


def test_cartesian_cutting_plane_fixed_z(spherical_ds):
    ds = spherical_ds
    normal = np.array([0.0, 0.0, 1.0])
    center = np.array([0.0, 0.0, 0.5])
    slc = ds.cartesian_cutting(normal, center)
    zvals = slc["index", "z_val"].to("code_length").d
    assert np.allclose(zvals, ds.quan(0.5, "code_length").d, atol=0.05)


@pytest.mark.mpl_image_compare
def test_vertical_slice_at_sphere_edge(spherical_ds):
    ds = spherical_ds
    normal = np.array([0.0, 1.0, 0.0])
    center = np.array([0.0, 0.75, 0.0])
    slc = ds.cartesian_cutting(normal, center)
    frb = slc.to_frb(2.0, 50)
    vals = frb["index", "z_val"].to("code_length")
    vals[~frb.get_mask(("index", "z_val"))] = np.nan

    f, axs = plt.subplots(1)
    axs.imshow(vals, origin="lower", extent=frb.bounds)
    return f


def test_cartesian_cutting_plane_with_axis_ordering():
    # check that slicing works with any axis order
    shp = (32, 32, 32)
    axes = ["r", "theta", "phi"]
    bbox_ranges = {"r": [0.0, 1.0], "theta": [0, np.pi], "phi": [0, 2 * np.pi]}

    # set the attributes for the plane, including a north vector found
    # for an arbitrary point on the plane.
    normal = np.array([1.0, 1.0, 1.0])
    center = np.array([0.0, 0.0, 0.0])
    x, y = 1.0, 1.0
    z = -x * normal[0] - y * normal[1]
    north_pt = np.array([x, y, z])
    assert np.dot(normal, north_pt) == 0.0  # just to be sure...

    frb_vals = []
    for axis_order in itertools.permutations(axes):
        bbox = np.zeros((3, 2))
        for i, ax in enumerate(axis_order):
            bbox[i, :] = bbox_ranges[ax]
        ds = _get_spherical_uniform_grid(shp, bbox, tuple(axis_order))
        slc = ds.cartesian_cutting(normal, center, north_vector=north_pt)
        frb = slc.to_frb(2.0, 50)
        vals = frb["index", "z_val"].to("code_length")
        vals[~frb.get_mask(("index", "z_val"))] = np.nan
        frb_vals.append(vals.d)

    for frb_z in frb_vals[1:]:
        np.allclose(frb_z, frb_vals[0])

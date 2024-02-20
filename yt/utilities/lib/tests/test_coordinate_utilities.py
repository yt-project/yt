import numpy as np

from yt.utilities.lib.coordinate_utilities import (
    CartesianMixedCoordBBox,
    SphericalMixedCoordBBox,
    cartesian_bboxes,
    cartesian_points_to_spherical,
    spherical_points_to_cartesian,
)


def test_cartesian_bboxes_for_spherical():

    # this test checks the special cases where
    # a spherical volume element crosses an axis
    # or when an element edge is lined up with an axis

    # check element that includes theta=0 as an edge
    r = np.array([0.95])
    dr = np.array([0.1])
    theta = np.array([0.05])
    dtheta = np.array([0.1])
    phi = np.array([0.05])
    dphi = np.array([0.05])

    x = np.full(r.shape, np.nan, dtype="float64")
    y = np.full(r.shape, np.nan, dtype="float64")
    z = np.full(r.shape, np.nan, dtype="float64")
    dx = np.full(r.shape, np.nan, dtype="float64")
    dy = np.full(r.shape, np.nan, dtype="float64")
    dz = np.full(r.shape, np.nan, dtype="float64")

    bbox_handler = SphericalMixedCoordBBox()

    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)
    assert z + dz / 2 == 1.0
    assert np.allclose(x - dx / 2, 0.0)
    assert np.allclose(y - dy / 2, 0.0)

    # now theta = np.pi
    theta = np.array([np.pi - dtheta[0] / 2])
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)
    assert z - dz / 2 == -1.0
    assert np.allclose(x - dx / 2, 0.0)
    assert np.allclose(y - dy / 2, 0.0)

    # element at equator, overlapping the +y axis
    theta = np.array([np.pi / 2])
    phi = np.array([np.pi / 2])
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)

    assert y + dy / 2 == 1.0
    assert np.allclose(x, 0.0)
    assert np.allclose(z, 0.0)

    # element at equator, overlapping the -x axis
    phi = np.array([np.pi])
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)

    assert x - dx / 2 == -1.0
    assert np.allclose(y, 0.0)
    assert np.allclose(z, 0.0)

    # element at equator, overlapping the -y axis
    phi = np.array([3 * np.pi / 2])
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)

    assert y - dy / 2 == -1.0
    assert np.allclose(x, 0.0)
    assert np.allclose(z, 0.0)

    # element at equator, overlapping +x axis
    phi = dphi / 2.0
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)
    assert x + dx / 2 == 1.0

    # element with edge on +x axis in -theta direction
    theta = np.pi / 2 - dtheta / 2
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)
    assert x + dx / 2 == 1.0

    # element with edge on +x axis in +theta direction
    theta = np.pi / 2 + dtheta / 2
    cartesian_bboxes(bbox_handler, r, theta, phi, dr, dtheta, dphi, x, y, z, dx, dy, dz)
    assert x + dx / 2 == 1.0

    # finally, check that things work OK with a wide range of
    # angles

    r_edges = np.linspace(0.4, 1.0, 10, dtype="float64")
    theta_edges = np.linspace(0, np.pi, 10, dtype="float64")
    phi_edges = np.linspace(0.0, 2 * np.pi, 10, dtype="float64")

    r = (r_edges[0:-1] + r_edges[1:]) / 2.0
    theta = (theta_edges[0:-1] + theta_edges[1:]) / 2.0
    phi = (phi_edges[0:-1] + phi_edges[1:]) / 2.0

    dr = r_edges[1:] - r_edges[:-1]
    dtheta = theta_edges[1:] - theta_edges[:-1]
    dphi = phi_edges[1:] - phi_edges[:-1]

    r_th_ph = np.meshgrid(r, theta, phi)
    d_r_th_ph = np.meshgrid(dr, dtheta, dphi)
    r_th_ph = [r_th_ph[i].ravel() for i in range(3)]
    d_r_th_ph = [d_r_th_ph[i].ravel() for i in range(3)]

    x_y_z = [np.full(r_th_ph[0].shape, np.nan, dtype="float64") for _ in range(3)]
    d_x_y_z = [np.full(r_th_ph[0].shape, np.nan, dtype="float64") for _ in range(3)]

    cartesian_bboxes(
        bbox_handler,
        r_th_ph[0],
        r_th_ph[1],
        r_th_ph[2],
        d_r_th_ph[0],
        d_r_th_ph[1],
        d_r_th_ph[2],
        x_y_z[0],
        x_y_z[1],
        x_y_z[2],
        d_x_y_z[0],
        d_x_y_z[1],
        d_x_y_z[2],
    )

    assert np.all(np.isfinite(x_y_z))
    assert np.all(np.isfinite(d_x_y_z))

    # and check the extents again for completeness...
    for i in range(3):
        max_val = np.max(x_y_z[i] + d_x_y_z[i] / 2.0)
        min_val = np.min(x_y_z[i] - d_x_y_z[i] / 2.0)
        assert max_val == 1.0
        assert min_val == -1.0


def test_cartesian_passthrough():

    bbox_handler = CartesianMixedCoordBBox()
    rng = np.random.default_rng()
    sz = 10

    xyz_in = [rng.random(sz) for _ in range(3)]
    dxyz_in = [rng.random(sz) for _ in range(3)]
    xyz_out = [np.full(xyz_in[0].shape, np.nan) for _ in range(3)]
    dxyz_out = [np.full(xyz_in[0].shape, np.nan) for _ in range(3)]

    cartesian_bboxes(
        bbox_handler,
        xyz_in[0],
        xyz_in[1],
        xyz_in[2],
        dxyz_in[0],
        dxyz_in[1],
        dxyz_in[2],
        xyz_out[0],
        xyz_out[1],
        xyz_out[2],
        dxyz_out[0],
        dxyz_out[1],
        dxyz_out[2],
    )

    assert np.all(np.isfinite(xyz_out))
    for idim in range(3):
        assert np.all(xyz_in[idim] == xyz_out[idim])
        assert np.all(dxyz_in[idim] == dxyz_out[idim])


def test_spherical_cartesian_roundtrip():
    xyz = [np.linspace(0, 1, 10) for _ in range(3)]
    xyz = np.meshgrid(*xyz)
    xyz = [xyzi.ravel() for xyzi in xyz]
    x, y, z = xyz

    r, theta, phi = cartesian_points_to_spherical(x, y, z)
    x_out, y_out, z_out = spherical_points_to_cartesian(r, theta, phi)

    assert np.allclose(x_out, x)
    assert np.allclose(y_out, y)
    assert np.allclose(z_out, z)
    assert np.max(r) == np.sqrt(3.0)

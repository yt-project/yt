import numpy as np

from yt.testing import assert_almost_equal, fake_sph_orientation_ds, requires_module
from yt.utilities.lib.pixelization_routines import pixelize_sph_kernel_projection
from yt.utilities.on_demand_imports import _scipy
from yt.visualization.volume_rendering import off_axis_projection as OffAP

spatial = _scipy.spatial
ndimage = _scipy.ndimage


def test_no_rotation():
    """Determines if a projection processed through
    off_axis_projection with no rotation will give the same
    image buffer if processed directly through
    pixelize_sph_kernel_projection
    """
    normal_vector = [0.0, 0.0, 1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    ad = ds.all_data()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = (left_edge + right_edge) / 2
    width = right_edge - left_edge
    px = ad[("all", "particle_position_x")]
    py = ad[("all", "particle_position_y")]
    hsml = ad[("all", "smoothing_length")]
    quantity_to_smooth = ad[("gas", "density")]
    density = ad[("io", "density")]
    mass = ad[("io", "particle_mass")]
    bounds = [-4, 4, -4, 4, -4, 4]

    buf2 = np.zeros(resolution)
    buf1 = OffAP.off_axis_projection(
        ds, center, normal_vector, width, resolution, ("gas", "density")
    )
    pixelize_sph_kernel_projection(
        buf2, px, py, hsml, mass, density, quantity_to_smooth, bounds
    )
    assert_almost_equal(buf1.ndarray_view(), buf2)


@requires_module("scipy")
def test_basic_rotation_1():
    """All particles on Z-axis should now be on the negative Y-Axis
    fake_sph_orientation has three z-axis particles,
    so there should be three y-axis particles after rotation
    (0, 0, 1) -> (0, -1)
    (0, 0, 2) -> (0, -2)
    (0, 0, 3) -> (0, -3)
    In addition, we should find a local maxima at (0, 0) due to:
    (0, 0, 0) -> (0, 0)
    (0, 1, 0) -> (0, 0)
    (0, 2, 0) -> (0, 0)
    and the one particle on the x-axis should not change its position:
    (1, 0, 0) -> (1, 0)
    """
    expected_maxima = ([0.0, 0.0, 0.0, 0.0, 1.0], [0.0, -1.0, -2.0, -3.0, 0.0])
    normal_vector = [0.0, 1.0, 0.0]
    north_vector = [0.0, 0.0, -1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = (left_edge + right_edge) / 2
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds,
        center,
        normal_vector,
        width,
        resolution,
        ("gas", "density"),
        north_vector=north_vector,
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_basic_rotation_2():
    """Rotation of x-axis onto z-axis.
    All particles on z-axis should now be on the negative x-Axis fake_sph_orientation
    has three z-axis particles, so there should be three x-axis particles after rotation
    (0, 0, 1) -> (-1, 0)
    (0, 0, 2) -> (-2, 0)
    (0, 0, 3) -> (-3, 0)
    In addition, we should find a local maxima at (0, 0) due to:
    (0, 0, 0) -> (0, 0)
    (1, 0, 0) -> (0, 0)
    and the two particles on the y-axis should not change its position:
    (0, 1, 0) -> (0, 1)
    (0, 2, 0) -> (0, 2)
    """
    expected_maxima = (
        [-1.0, -2.0, -3.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, 2.0],
    )
    normal_vector = [1.0, 0.0, 0.0]
    north_vector = [0.0, 1.0, 0.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = (left_edge + right_edge) / 2
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds,
        center,
        normal_vector,
        width,
        resolution,
        ("gas", "density"),
        north_vector=north_vector,
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_basic_rotation_3():
    """Rotation of z-axis onto negative z-axis.
    All fake particles on z-axis should now be of the negative z-Axis.
    fake_sph_orientation has three z-axis particles,
    so we should have a local maxima at (0, 0)
    (0, 0, 1) -> (0, 0)
    (0, 0, 2) -> (0, 0)
    (0, 0, 3) -> (0, 0)
    In addition, (0, 0, 0) should also contribute to the local maxima at (0, 0):
    (0, 0, 0) -> (0, 0)
    x-axis particles should be rotated as such:
    (1, 0, 0) -> (0, -1)
    and same goes for y-axis particles:
    (0, 1, 0) -> (-1, 0)
    (0, 2, 0) -> (-2, 0)
    """
    expected_maxima = ([0.0, 0.0, -1.0, -2.0], [0.0, -1.0, 0.0, 0.0])
    normal_vector = [0.0, 0.0, -1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = (left_edge + right_edge) / 2
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds, center, normal_vector, width, resolution, ("gas", "density")
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_basic_rotation_4():
    """Rotation of x-axis to z-axis and original z-axis to y-axis with the use
    of the north_vector. All fake particles on z-axis should now be on the
    y-Axis.  All fake particles on the x-axis should now be on the z-axis, and
    all fake particles on the y-axis should now be on the x-axis.

    (0, 0, 1) -> (0, 1)
    (0, 0, 2) -> (0, 2)
    (0, 0, 3) -> (0, 3)
    In addition, (0, 0, 0) should contribute to the local maxima at (0, 0):
    (0, 0, 0) -> (0, 0)
    x-axis particles should be rotated and contribute to the local maxima at (0, 0):
    (1, 0, 0) -> (0, 0)
    and the y-axis particles shift into the positive x direction:
    (0, 1, 0) -> (1, 0)
    (0, 2, 0) -> (2, 0)
    """
    expected_maxima = ([0.0, 0.0, 0.0, 0.0, 1.0, 2.0], [1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
    normal_vector = [1.0, 0.0, 0.0]
    north_vector = [0.0, 0.0, 1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = (left_edge + right_edge) / 2
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds,
        center,
        normal_vector,
        width,
        resolution,
        ("gas", "density"),
        north_vector=north_vector,
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_center_1():
    """Change the center to [0, 3, 0]
    Every point will be shifted by 3 in the y-domain
    With this, we should not be able to see any of the y-axis particles
    (0, 1, 0) -> (0, -2)
    (0, 2, 0) -> (0, -1)
    (0, 0, 1) -> (0, -3)
    (0, 0, 2) -> (0, -3)
    (0, 0, 3) -> (0, -3)
    (0, 0, 0) -> (0, -3)
    (1, 0, 0) -> (1, -3)
    """
    expected_maxima = ([0.0, 0.0, 0.0, 1.0], [-2.0, -1.0, -3.0, -3.0])
    normal_vector = [0.0, 0.0, 1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    # center = [(left_edge[0] + right_edge[0])/2,
    #            left_edge[1],
    #           (left_edge[2] + right_edge[2])/2]
    center = [0.0, 3.0, 0.0]
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds, center, normal_vector, width, resolution, ("gas", "density")
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_center_2():
    """Change the center to [0, -1, 0]
    Every point will be shifted by 1 in the y-domain
    With this, we should not be able to see any of the y-axis particles
    (0, 1, 0) -> (0, 2)
    (0, 2, 0) -> (0, 3)
    (0, 0, 1) -> (0, 1)
    (0, 0, 2) -> (0, 1)
    (0, 0, 3) -> (0, 1)
    (0, 0, 0) -> (0, 1)
    (1, 0, 0) -> (1, 1)
    """
    expected_maxima = ([0.0, 0.0, 0.0, 1.0], [2.0, 3.0, 1.0, 1.0])
    normal_vector = [0.0, 0.0, 1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [0.0, -1.0, 0.0]
    width = right_edge - left_edge
    buf1 = OffAP.off_axis_projection(
        ds, center, normal_vector, width, resolution, ("gas", "density")
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def test_center_3():
    """Change the center to the left edge, or [0, -8, 0]
    Every point will be shifted by 8 in the y-domain
    With this, we should not be able to see anything !
    """
    expected_maxima = ([], [])
    normal_vector = [0.0, 0.0, 1.0]
    resolution = (64, 64)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [0.0, -1.0, 0.0]
    width = [
        (right_edge[0] - left_edge[0]),
        left_edge[1],
        (right_edge[2] - left_edge[2]),
    ]
    buf1 = OffAP.off_axis_projection(
        ds, center, normal_vector, width, resolution, ("gas", "density")
    )
    find_compare_maxima(expected_maxima, buf1, resolution, width)


@requires_module("scipy")
def find_compare_maxima(expected_maxima, buf, resolution, width):
    buf_ndarray = buf.ndarray_view()
    max_filter_buf = ndimage.maximum_filter(buf_ndarray, size=5)
    maxima = np.isclose(max_filter_buf, buf_ndarray, rtol=1e-09)

    # ignore contributions from zones of no smoothing
    for i in range(len(maxima)):
        for j in range(len(maxima[i])):
            if np.isclose(buf_ndarray[i, j], 0.0, 1e-09):
                maxima[i, j] = False
    coords = ([], [])

    for i in range(len(maxima)):
        for j in range(len(maxima[i])):
            if maxima[i, j]:
                coords[0].append(i)
                coords[1].append(j)
    pixel_tolerance = 2.0
    x_scaling_factor = resolution[0] / width[0]
    y_scaling_factor = resolution[1] / width[1]
    for i in range(len(expected_maxima[0])):
        found_match = False
        for j in range(len(coords[0])):
            # normalize coordinates
            x_coord = coords[0][j]
            y_coord = coords[1][j]
            x_coord -= resolution[0] / 2
            y_coord -= resolution[1] / 2
            x_coord /= x_scaling_factor
            y_coord /= y_scaling_factor
            if (
                spatial.distance.euclidean(
                    [x_coord, y_coord], [expected_maxima[0][i], expected_maxima[1][i]]
                )
                < pixel_tolerance
            ):
                found_match = True
                break
        if found_match is not True:
            raise AssertionError
    pass

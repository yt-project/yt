import yt
from yt.visualization.volume_rendering import \
    off_axis_projection as OffAP
from yt.testing import \
    fake_sph_orientation_ds, \
    assert_equal, \
    assert_almost_equal, \
    requires_module
from yt.utilities.lib.pixelization_routines import \
    pixelize_sph_kernel_projection
import numpy as np
from scipy.ndimage.filters import \
    maximum_filter, minimum_filter
from scipy.spatial.distance import\
    euclidean
import matplotlib.pyplot as plt

def test_no_rotation():
    """ Determines if a projection processed through 
    off_axis_projection with no rotation will give the same
    image buffer if processed directly through 
    pixelize_sph_kernel_projection
    """

    normal_vector = [0., 0., 1.]
    resolution = (128, 128)
    ds = fake_sph_orientation_ds()
    ad = ds.all_data()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [(left_edge[0] + right_edge[0])/2,
              (left_edge[1] + right_edge[1])/2,
              (left_edge[2] + right_edge[2])/2]
    width = [(right_edge[0] - left_edge[0]),
             (right_edge[1] - left_edge[1]),
             (right_edge[2] - left_edge[2])]
    px = ad["particle_position_x"]
    py = ad["particle_position_y"]
    hsml = ad["smoothing_length"]
    quantity_to_smooth = ad["particle_mass"]
    density = ad["density"]
    mass = ad["particle_mass"]
    bounds = [-4, 4, -4, 4, -4, 4]

    buf2 = np.zeros(resolution)
    buf1 = OffAP.off_axis_projection(ds,
                                     center,
                                     normal_vector,
                                     width,
                                     resolution,
                                     'particle_mass'
                                     )
    pixelize_sph_kernel_projection(buf2,
                                   px,
                                   py,
                                   hsml,
                                   mass,
                                   density,
                                   quantity_to_smooth,
                                   bounds
                                   )
    assert_almost_equal(buf1, buf2)

@requires_module('scipy')
def test_basic_rotation_1():
    """ Rotation of z-axis onto y-axis. All fake particles on Z-axis should now be on the Y-Axis
        fake_sph_orientation has three z-axis particles, so there should be three y-axis particles 
        after rotation 
        (0, 0, 1) -> (0, 1)
        (0, 0, 2) -> (0, 2)
        (0, 0, 3) -> (0, 3)
        In addition, we should find a local maxima at (0, 0) due to:
        (0, 0, 0) -> (0, 0)
        (0, 1, 0) -> (0, 0)
        (0, 2, 0) -> (0, 0)
        and the one particle on the x-axis should not change its position:
        (1, 0, 0) -> (1, 0)
    """
    expected_maxima = ([0., 0., 0., 0., 1.], [0., 1., 2., 3., 0.])
    normal_vector = [0., 1., 0.]
    resolution = (256, 256)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [(left_edge[0] + right_edge[0])/2,
              (left_edge[1] + right_edge[1])/2,
              (left_edge[2] + right_edge[2])/2]
    width = [(right_edge[0] - left_edge[0]),
             (right_edge[1] - left_edge[1]),
             (right_edge[2] - left_edge[2])]    
    buf1 = OffAP.off_axis_projection(ds,
                                     center,
                                     normal_vector,
                                     width,
                                     resolution,
                                     'particle_mass'
                                     )
    max_filter_buf = maximum_filter(buf1, size=5)
    maxima = np.isclose(max_filter_buf, buf1, rtol=1e-09)
    # plt.imsave('sph_image.png', buf1)
    # plt.imsave('max_filter.png', max_filter_buf)
    # ignore contributions from zones of no smoothing
    for i in range(len(maxima)):
        for j in range(len(maxima[i])):
            if np.isclose(buf1[i, j], 0., 1e-09):
                maxima[i, j] = False
    coords = ([], [])
    # Using a step size two since the same maxima is often double counted
    for i in range(0, len(maxima), 2):
        for j in range(0, len(maxima[i]), 2):
            if maxima[i, j] == True:
                coords[0].append(i)
                coords[1].append(j)
    if len(coords[0]) == 0:
        assert False
    #assert that our expected maxima coordinates are in fact in the image
    pixel_tolerance = 2.0
    x_scaling_factor = resolution[0] / width[0]
    y_scaling_factor = resolution[1] / width[1]
    for i in range(len(coords[0])):
        #normalize coordinates
        x_coord = coords[0][i]
        y_coord = coords[1][i]
        x_coord -= resolution[0] / 2
        y_coord -= resolution[0] / 2
        found_match = False
        for j in range(len(expected_maxima[0])):
            if euclidean([x_coord, y_coord], [x_scaling_factor * expected_maxima[0][j], 
                                              y_scaling_factor * expected_maxima[1][j]]) < pixel_tolerance:
                found_match = True 
                break
        if found_match is not True:
            assert False
    pass

@requires_module('scipy')
def test_basic_rotation_2():
    """ Rotation of z-axis onto x-axis. All fake particles on z-axis should now be on the x-Axis
    fake_sph_orientation has three z-axis particles, so there should be three x-axis particles 
    after rotation 
    (0, 0, 1) -> (1, 0)
    (0, 0, 2) -> (2, 0)
    (0, 0, 3) -> (3, 0)
    In addition, we should find a local maxima at (0, 0) due to:
    (0, 0, 0) -> (0, 0)
    (1, 0, 0) -> (0, 0)
    and the two particles on the y-axis should not change its position:
    (0, 1, 0) -> (0, 1)
    (0, 2, 0) -> (0, 2)
    """ 
    expected_maxima = ([1., 2., 3., 0., 0., 0.], 
                       [0., 0., 0., 0., 1., 2.])
    normal_vector = [1., 0., 0.]
    resolution = (256, 256)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [(left_edge[0] + right_edge[0])/2,
              (left_edge[1] + right_edge[1])/2,
              (left_edge[2] + right_edge[2])/2]
    width = [(right_edge[0] - left_edge[0]),
             (right_edge[1] - left_edge[1]),
             (right_edge[2] - left_edge[2])]    
    buf1 = OffAP.off_axis_projection(ds,
                                     center,
                                     normal_vector,
                                     width,
                                     resolution,
                                     'particle_mass'
                                     )
    max_filter_buf = maximum_filter(buf1, size=5)
    maxima = np.isclose(max_filter_buf, buf1, rtol=1e-09)
    # ignore contributions from zones of no smoothing
    for i in range(len(maxima)):
        for j in range(len(maxima[i])):
            if np.isclose(buf1[i, j], 0., 1e-09):
                maxima[i, j] = False
    coords = ([], [])

    # Using a step size two since the same maxima is often double counted
    for i in range(0, len(maxima), 2):
        for j in range(0, len(maxima[i]), 2):
            if maxima[i, j] == True:
                coords[0].append(i)
                coords[1].append(j)
    #assert that our expected maxima coordinates are in fact in the image
    pixel_tolerance = 2.0
    x_scaling_factor = resolution[0] / width[0]
    y_scaling_factor = resolution[1] / width[1]
    for i in range(len(coords[0])):
        #normalize coordinates
        x_coord = coords[0][i]
        y_coord = coords[1][i]
        x_coord -= resolution[0] / 2
        y_coord -= resolution[0] / 2
        found_match = False
        for j in range(len(expected_maxima[0])):
            if euclidean([x_coord, y_coord], [x_scaling_factor * expected_maxima[0][j], 
                                              y_scaling_factor * expected_maxima[1][j]]) < pixel_tolerance:
                found_match = True 
                break
        if found_match is not True:
            assert False
    pass


@requires_module('scipy')
def test_basic_rotation_3():
    """ Rotation of z-axis onto negative z-axis. All fake particles on z-axis should now be on 
    the negative z-Axis.
    fake_sph_orientation has three z-axis particles, so we should have a local maxima at (0, 0)
    (0, 0, 1) -> (0, 0)
    (0, 0, 2) -> (0, 0)
    (0, 0, 3) -> (0, 0)
    In addition, (0, 0, 0) should also contribute to the local maxima at (0, 0):
    (0, 0, 0) -> (0, 0)
    x-axis particles should be rotated 180 degrees
    (1, 0, 0) -> (0, -1)
    and same goes for y-axis particles:
    (0, 1, 0) -> (-1, 0)
    (0, 2, 0) -> (-2, 0)
    """ 
    expected_maxima = ([0., 0., -1., -2.], [0., -1., 0., 0.])
    normal_vector = [0., 0., -1.]
    resolution = (256, 256)
    ds = fake_sph_orientation_ds()
    left_edge = ds.domain_left_edge
    right_edge = ds.domain_right_edge
    center = [(left_edge[0] + right_edge[0])/2,
              (left_edge[1] + right_edge[1])/2,
              (left_edge[2] + right_edge[2])/2]
    width = [(right_edge[0] - left_edge[0]),
             (right_edge[1] - left_edge[1]),
             (right_edge[2] - left_edge[2])]    
    buf1 = OffAP.off_axis_projection(ds,
                                     center,
                                     normal_vector,
                                     width,
                                     resolution,
                                     'particle_mass'
                                     )
    max_filter_buf = maximum_filter(buf1, size=5)
    maxima = np.isclose(max_filter_buf, buf1, rtol=1e-09)
    plt.imsave('sph_image.png', buf1)
    plt.imsave('max_filter.png', max_filter_buf)
    # ignore contributions from zones of no smoothing
    for i in range(len(maxima)):
        for j in range(len(maxima[i])):
            if np.isclose(buf1[i, j], 0., 1e-09):
                maxima[i, j] = False
    coords = ([], [])
    # Using a step size two since the same maxima is often double/quadruple counted
    for i in range(0, len(maxima), 2):
        for j in range(0, len(maxima[i]), 2):
            if maxima[i, j] == True:
                coords[0].append(i)
                coords[1].append(j)

    #assert that our expected maxima coordinates are in fact in the image
    pixel_tolerance = 2.0
    x_scaling_factor = resolution[0] / width[0]
    y_scaling_factor = resolution[1] / width[1]
    for i in range(len(coords[0])):
        #normalize coordinates
        x_coord = coords[0][i]
        y_coord = coords[1][i]
        x_coord -= resolution[0] / 2
        y_coord -= resolution[0] / 2
        found_match = False
        for j in range(len(expected_maxima[0])):
            if euclidean([x_coord, y_coord], [x_scaling_factor * expected_maxima[0][j], 
                                              y_scaling_factor * expected_maxima[1][j]]) < pixel_tolerance:
                found_match = True 
                break
        if found_match is not True:
            assert False
    pass
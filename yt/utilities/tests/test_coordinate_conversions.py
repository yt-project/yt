import numpy as np

from yt.testing import assert_array_almost_equal
from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component, \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, get_sph_r, \
    get_sph_theta, get_sph_phi

# Randomly generated coordinates in the domain [[-1,1],[-1,1],-1,1]]
coords = np.array([[-0.41503037, -0.22102472, -0.55774212],
                   [ 0.73828247, -0.17913899,  0.64076921],
                   [ 0.08922066, -0.94254844, -0.61774511],
                   [ 0.10173242, -0.95789145,  0.16294352],
                   [ 0.73186508, -0.3109153 ,  0.75728738],
                   [ 0.8757989 , -0.41475119, -0.57039201],
                   [ 0.58040762,  0.81969082,  0.46759728],
                   [-0.89983356, -0.9853683 , -0.38355343]]).T

def test_spherical_coordinate_conversion():
    normal = [0, 0, 1]
    real_r =     [ 0.72950559,  0.99384957,  1.13047198,  0.97696269,  
                   1.09807968,  1.12445067,  1.10788685,  1.38843954]
    real_theta = [ 2.44113629,  0.87012028,  2.14891444,  1.4032274 ,  
                   0.80979483,  2.10280198,  1.13507735,  1.85068416]
    real_phi =   [-2.65224483, -0.23804243, -1.47641858, -1.46498842, 
                  -0.40172325, -0.4422801 ,  0.95466734, -2.31085392]

    calc_r = get_sph_r(coords)
    calc_theta = get_sph_theta(coords, normal)
    calc_phi = get_sph_phi(coords, normal)

    assert_array_almost_equal(calc_r, real_r)
    assert_array_almost_equal(calc_theta, real_theta)
    assert_array_almost_equal(calc_phi, real_phi)

    normal = [1, 0, 0]
    real_theta = [ 2.17598842,  0.73347681,  1.49179079,  1.46647589,  
                   0.8412984 ,  0.67793705,  1.0193883 ,  2.27586987]
    real_phi =   [-0.37729951, -2.86898397, -0.99063518, -1.73928995, 
                   -2.75201227,-0.62870527,  2.08920872, -1.19959244]

    calc_theta = get_sph_theta(coords, normal)
    calc_phi = get_sph_phi(coords, normal)
    
    assert_array_almost_equal(calc_theta, real_theta)
    assert_array_almost_equal(calc_phi, real_phi)

def test_cylindrical_coordiante_conversion():
    normal = [0, 0, 1]
    real_r =     [ 0.47021498,  0.75970506,  0.94676179,  0.96327853,  
                   0.79516968,  0.96904193,  1.00437346,  1.3344104 ]    
    real_theta = [-2.65224483, -0.23804243, -1.47641858, -1.46498842, 
                  -0.40172325, -0.4422801 ,  0.95466734, -2.31085392]
    real_z =     [-0.55774212,  0.64076921, -0.61774511,  0.16294352,
                   0.75728738, -0.57039201,  0.46759728, -0.38355343]

    calc_r = get_cyl_r(coords, normal)
    calc_theta = get_cyl_theta(coords, normal)
    calc_z = get_cyl_z(coords, normal)

    assert_array_almost_equal(calc_r, real_r)
    assert_array_almost_equal(calc_theta, real_theta)
    assert_array_almost_equal(calc_z, real_z)

    normal = [1, 0, 0]
    real_r =     [ 0.59994016,  0.66533898,  1.12694569,  0.97165149,
                   0.81862843,  0.70524152,  0.94368441,  1.05738542]
    real_theta = [-0.37729951, -2.86898397, -0.99063518, -1.73928995, 
                  -2.75201227, -0.62870527,  2.08920872, -1.19959244]
    real_z =     [-0.41503037,  0.73828247,  0.08922066,  0.10173242,
                   0.73186508,  0.8757989 ,  0.58040762, -0.89983356]

    calc_r = get_cyl_r(coords, normal)
    calc_theta = get_cyl_theta(coords, normal)
    calc_z = get_cyl_z(coords, normal)

    assert_array_almost_equal(calc_r, real_r)
    assert_array_almost_equal(calc_theta, real_theta)
    assert_array_almost_equal(calc_z, real_z)

def test_spherical_coordinate_projections():
    normal = [0, 0, 1]
    theta = get_sph_theta(coords, normal)
    phi = get_sph_phi(coords, normal)
    zero = np.tile(0,coords.shape[1])

    # Purely radial field
    vecs = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    assert_array_almost_equal(zero, get_sph_theta_component(vecs, theta, phi, normal))
    assert_array_almost_equal(zero, get_sph_phi_component(vecs, phi, normal))

    # Purely toroidal field
    vecs = np.array([-np.sin(phi), np.cos(phi), zero])
    assert_array_almost_equal(zero, get_sph_theta_component(vecs, theta, phi, normal))
    assert_array_almost_equal(zero, get_sph_r_component(vecs, theta, phi, normal))

    # Purely poloidal field
    vecs = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)])
    assert_array_almost_equal(zero, get_sph_phi_component(vecs, phi, normal))
    assert_array_almost_equal(zero, get_sph_r_component(vecs, theta, phi, normal))

def test_cylindrical_coordinate_projections():
    normal = [0, 0, 1]
    theta = get_cyl_theta(coords, normal)
    z = get_cyl_z(coords, normal)
    zero = np.tile(0, coords.shape[1])

    # Purely radial field
    vecs = np.array([np.cos(theta), np.sin(theta), zero])
    assert_array_almost_equal(zero, get_cyl_theta_component(vecs, theta, normal))
    assert_array_almost_equal(zero, get_cyl_z_component(vecs, normal))

    # Purely toroidal field
    vecs = np.array([-np.sin(theta), np.cos(theta), zero])
    assert_array_almost_equal(zero, get_cyl_z_component(vecs, normal))
    assert_array_almost_equal(zero, get_cyl_r_component(vecs, theta, normal))

    # Purely z field
    vecs = np.array([zero, zero, z])
    assert_array_almost_equal(zero, get_cyl_theta_component(vecs, theta, normal))
    assert_array_almost_equal(zero, get_cyl_r_component(vecs, theta, normal))

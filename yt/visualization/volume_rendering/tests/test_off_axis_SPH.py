import yt
from yt.visualization.volume_rendering import \
    off_axis_projection as OffAP
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def test_create_projection():
    print('running test')
    ds = yt.load("../../../../../Data/GadgetDiskGalaxy/snapshot_200.hdf5")
    # ds = yt.load('../../../../../Data/IsothermalCollapse/snap_505', unit_base={'UnitLength_in_cm': 5.0e16, 
    #              'UnitMass_in_g': 1.98992e33, 'UnitVelocity_in_cm_per_s': 46385.190}, 
    #              bounding_box = [[-3, 3], [-3, 3], [-3, 3]])    
    resolution = (512, 512)
    normal_vector = np.array([0., 0., 1.])
    buf = OffAP.off_axis_projection(ds,
                              [32000., 32000., 32000.],
                              normal_vector,
                              [64000., 64000., 64000.],
                              resolution,
                              'particle_mass')
    np.seterr(divide='ignore')
    #plt.imsave('sph_img.png', np.log10(buf))
    plt.imsave('sph_img.png', buf)

def make_gif():
    normal_vector = np.array([-2., 2., -5])
    ds = yt.load('../../../../../Data/IsothermalCollapse/snap_505', unit_base={'UnitLength_in_cm': 5.0e16, 
                 'UnitMass_in_g': 1.98992e33, 'UnitVelocity_in_cm_per_s': 46385.190}, 
                 bounding_box = [[-3, 3], [-3, 3], [-3, 3]])    
    resolution = (512, 512)
    for i in range(160):
        normal_vector[2] += 0.0625
        buf = OffAP.off_axis_projection(ds,
                              [0., 0., 0.],
                              normal_vector,
                              [6., 6., 6.],
                              resolution,
                              'particle_mass')
        plt.imsave('images/img_' + str(i) + '.png',
                   buf)
        print('img_' + str(i) + '.png')
if __name__ == "__main__":
    test_create_projection()
    #make_gif()

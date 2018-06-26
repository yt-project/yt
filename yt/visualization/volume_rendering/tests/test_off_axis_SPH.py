import yt
from yt.visualization.volume_rendering import \
    off_axis_projection as OffAP
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def test_create_projection():
    print('running test')
    ds = yt.load("../../../../../Data/GadgetDiskGalaxy/snapshot_200.hdf5")
    # ds = yt.load('../Data/IsothermalCollapse/snap_505', unit_base={'UnitLength_in_cm': 5.0e16, 
    #              'UnitMass_in_g': 1.98992e33, 'UnitVelocity_in_cm_per_s': 46385.190}, 
    #              bounding_box = [[-3, 3], [-3, 3], [-3, 3]])    
    resolution = (512, 512)
    normal_vector = np.array([0., 0., 1.])
    buf = OffAP.off_axis_projection(ds,
                              [0, 0, 0],
                              normal_vector,
                              0,
                              resolution,
                              'particle_mass')
    np.seterr(divide='ignore')
    #plt.imsave('sph_img.png', np.log10(buf))
    plt.imsave('sph_img.png', buf)

    
if __name__ == "__main__":
    test_create_projection()

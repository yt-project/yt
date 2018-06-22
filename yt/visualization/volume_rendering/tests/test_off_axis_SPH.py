import yt
from yt.visualization.volume_rendering import \
    off_axis_projection as OffAP
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def test_create_projection():
    print('running test')
    ds = yt.load("../Data/GadgetDiskGalaxy/snapshot_200.hdf5")
    #ad = ds.all_data()
    resolution = (512, 512)
    normal_vector = [0, 0, 1]
    buf = OffAP.off_axis_projection(ds,
                              [0, 0, 0],
                              normal_vector,
                              0,
                              resolution,
                              'particle_mass')
    np.seterr(divide='ignore')
    plt.imsave('sph_img.png', np.log10(buf))
    #plt.imsave('sph_img.png', buf)

    
if __name__ == "__main__":
    test_create_projection()

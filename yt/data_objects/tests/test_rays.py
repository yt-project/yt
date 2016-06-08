import numpy as np

from yt import \
    load
from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    assert_rel_equal, \
    requires_file
from yt.units.yt_array import \
    uconcatenate

def test_ray():
    for nproc in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs=nproc)
        dx = (ds.domain_right_edge - ds.domain_left_edge) / \
          ds.domain_dimensions
        # Three we choose, to get varying vectors, and ten random
        pp1 = np.random.random((3, 13))
        pp2 = np.random.random((3, 13))
        pp1[:,0] = [0.1, 0.2, 0.3]
        pp2[:,0] = [0.8, 0.1, 0.4]
        pp1[:,1] = [0.9, 0.2, 0.3]
        pp2[:,1] = [0.8, 0.1, 0.4]
        pp1[:,2] = [0.9, 0.2, 0.9]
        pp2[:,2] = [0.8, 0.1, 0.4]
        unitary = ds.arr(1.0, '')
        for i in range(pp1.shape[1]):
            p1 = ds.arr(pp1[:,i] + 1e-8 * np.random.random(3), 'code_length')
            p2 = ds.arr(pp2[:,i] + 1e-8 * np.random.random(3), 'code_length')

            my_ray = ds.ray(p1, p2)
            yield assert_rel_equal, my_ray['dts'].sum(), unitary, 14
            ray_cells = my_ray['dts'] > 0

            # find cells intersected by the ray
            my_all = ds.all_data()
            
            dt = np.abs(dx / (p2 - p1))
            tin  = uconcatenate([[(my_all['x'] - p1[0]) / (p2 - p1)[0] - 0.5 * dt[0]],
                                 [(my_all['y'] - p1[1]) / (p2 - p1)[1] - 0.5 * dt[1]],
                                 [(my_all['z'] - p1[2]) / (p2 - p1)[2] - 0.5 * dt[2]]])
            tout = uconcatenate([[(my_all['x'] - p1[0]) / (p2 - p1)[0] + 0.5 * dt[0]],
                                 [(my_all['y'] - p1[1]) / (p2 - p1)[1] + 0.5 * dt[1]],
                                 [(my_all['z'] - p1[2]) / (p2 - p1)[2] + 0.5 * dt[2]]])
            tin = tin.max(axis=0)
            tout = tout.min(axis=0)
            my_cells = (tin < tout) & (tin < 1) & (tout > 0)

            yield assert_equal, ray_cells.sum(), my_cells.sum()
            yield assert_rel_equal, my_ray['density'][ray_cells].sum(), \
                                    my_all['density'][my_cells].sum(), 14
            yield assert_rel_equal, my_ray['dts'].sum(), unitary, 14

@requires_file('GadgetDiskGalaxy/snapshot_200.hdf5')
def test_ray_in_particle_octree():
    ds = load('GadgetDiskGalaxy/snapshot_200.hdf5')
    start = ds.arr([31995.63476562, 31473.6640625, 28969.88671875], "code_length")
    end = ds.arr([31995.63476562, 31473.6640625, 29219.88671875], "code_length")
    ray = ds.ray(start, end)
    ray["t"]
    assert_equal(ray["dts"].sum(dtype="f8"), 1.0)

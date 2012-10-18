from yt.testing import *

def test_ray():
    pf = fake_random_pf(64, nprocs=8)
    dx = (pf.domain_right_edge - pf.domain_left_edge) / \
      pf.domain_dimensions

    p1 = np.random.random(3)
    p2 = np.random.random(3)

    my_ray = pf.h.ray(p1, p2)
    assert_equal(my_ray['dts'].sum(), 1.0)
    ray_cells = my_ray['dts'] > 0

    # find cells intersected by the ray
    my_all = pf.h.all_data()
    
    dt = np.abs(dx / (p2 - p1))
    tin  = np.concatenate([[(my_all['x'] - p1[0]) / (p2 - p1)[0] - 0.5 * dt[0]],
                           [(my_all['y'] - p1[1]) / (p2 - p1)[1] - 0.5 * dt[1]],
                           [(my_all['z'] - p1[2]) / (p2 - p1)[2] - 0.5 * dt[2]]])
    tout = np.concatenate([[(my_all['x'] - p1[0]) / (p2 - p1)[0] + 0.5 * dt[0]],
                           [(my_all['y'] - p1[1]) / (p2 - p1)[1] + 0.5 * dt[1]],
                           [(my_all['z'] - p1[2]) / (p2 - p1)[2] + 0.5 * dt[2]]])
    tin = tin.max(axis=0)
    tout = tout.min(axis=0)
    my_cells = (tin < tout) & (tin < 1) & (tout > 0)

    assert_equal(ray_cells.sum(), my_cells.sum())
    assert_equal(my_ray['Density'][ray_cells].sum(),
                 my_all['Density'][my_cells].sum())           

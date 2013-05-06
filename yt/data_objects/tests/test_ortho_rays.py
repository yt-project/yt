from yt.testing import *

def test_ortho_ray():
    pf = fake_random_pf(64, nprocs=8)
    dx = (pf.domain_right_edge - pf.domain_left_edge) / \
          pf.domain_dimensions

    axes = ['x', 'y', 'z']
    for ax, an in enumerate(axes):
        ocoord = np.random.random(2)

        my_oray = pf.h.ortho_ray(ax, ocoord)

        my_axes = range(3)
        del my_axes[ax]

        # find the cells intersected by the ortho ray
        my_all = pf.h.all_data()
        my_cells = (np.abs(my_all[axes[my_axes[0]]] - ocoord[0]) <= 
                    0.5 * dx[my_axes[0]]) & \
                   (np.abs(my_all[axes[my_axes[1]]] - ocoord[1]) <= 
                    0.5 * dx[my_axes[1]])

        yield assert_equal, my_oray['density'].sum(), \
                            my_all['density'][my_cells].sum()

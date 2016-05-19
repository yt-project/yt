import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_equal


def test_ortho_ray():
    ds = fake_random_ds(64, nprocs=8)
    dx = (ds.domain_right_edge - ds.domain_left_edge) / \
          ds.domain_dimensions

    axes = ['x', 'y', 'z']
    for ax, an in enumerate(axes):
        ocoord = ds.arr(np.random.random(2), 'code_length')

        my_oray = ds.ortho_ray(ax, ocoord)

        my_axes = ds.coordinates.x_axis[ax], ds.coordinates.y_axis[ax]

        # find the cells intersected by the ortho ray
        my_all = ds.all_data()
        my_cells = (np.abs(my_all[axes[my_axes[0]]] - ocoord[0]) <= 
                    0.5 * dx[my_axes[0]]) & \
                   (np.abs(my_all[axes[my_axes[1]]] - ocoord[1]) <= 
                    0.5 * dx[my_axes[1]])

        yield assert_equal, my_oray['density'].sum(), \
                            my_all['density'][my_cells].sum()

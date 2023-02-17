import numpy as np

import yt


def _construct_ds_from_coords_and_conn(coords, conn):
    # the distance from the origin
    node_data = {}
    dist = np.sum(coords**2, 1)
    node_data[("connect1", "test")] = dist[conn - 1]

    ds = yt.load_unstructured_mesh(conn - 1, coords, node_data=node_data)
    return ds


def test_mesh_slice_with_null_elems():
    # first construct a ds with all valid elements to get the count of
    # finite pixels in the final image
    coords = np.array(
        [
            [0.0, 1.0],
            [0.22746729, 0.0],
            [1.0, 0.87744062],
            [0.11373325, 0.5],
            [0.61373325, 0.43872031],
            [0.5, 0.93872031],
            [0.4, 0.6],
        ]
    )

    conn = np.array([[1, 2, 3, 4, 5, 6]])

    # the distance from the origin
    ds = _construct_ds_from_coords_and_conn(coords, conn)
    slc = yt.SlicePlot(ds, "z", ("connect1", "test"), buff_size=(10, 10))
    slc._setup_plots()
    expected_finite_pixels = np.isfinite(slc.frb[("connect1", "test")]).sum()

    # add on a null element, check that we get the same finite pixels
    null_element = [
        [
            7,
        ]
        * 6,
    ]
    conn = np.concatenate([conn, null_element])
    ds = _construct_ds_from_coords_and_conn(coords, conn)
    slc = yt.SlicePlot(ds, "z", ("connect1", "test"), buff_size=(10, 10))
    slc._setup_plots()
    actual_finite_pixels = np.isfinite(slc.frb[("connect1", "test")]).sum()

    assert expected_finite_pixels == actual_finite_pixels

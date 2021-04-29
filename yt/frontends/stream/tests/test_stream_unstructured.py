import numpy as np

from yt import SlicePlot, load_unstructured_mesh


def test_multi_mesh():
    coordsMulti = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=np.float64
    )

    connect1 = np.array(
        [
            [0, 1, 3],
        ],
        dtype=np.int64,
    )
    connect2 = np.array(
        [
            [1, 2, 3],
        ],
        dtype=np.int64,
    )

    data1 = {}
    data2 = {}
    data1["connect1", "test"] = np.array(
        [
            [0.0, 1.0, 3.0],
        ],
        dtype=np.float64,
    )
    data2["connect2", "test"] = np.array(
        [
            [1.0, 2.0, 3.0],
        ],
        dtype=np.float64,
    )

    connectList = [connect1, connect2]
    dataList = [data1, data2]

    ds = load_unstructured_mesh(connectList, coordsMulti, dataList)

    sl = SlicePlot(ds, "z", ("connect1", "test"))
    assert sl.data_source.field_data["connect1", "test"].shape == (1, 3)
    sl = SlicePlot(ds, "z", ("connect2", "test"))
    assert sl.data_source.field_data["connect2", "test"].shape == (1, 3)
    sl = SlicePlot(ds, "z", ("all", "test"))
    assert sl.data_source.field_data["all", "test"].shape == (2, 3)
    sl.annotate_mesh_lines()


def test_multi_field():
    coords = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], dtype=np.float64
    )

    connect = np.array([[0, 1, 3], [1, 2, 3]], dtype=np.int64)

    data = {}
    data["connect1", "test"] = np.array(
        [[0.0, 1.0, 3.0], [1.0, 2.0, 3.0]], dtype=np.float64
    )
    data["connect1", "testAgain"] = np.array(
        [[0.0, 1.0, 3.0], [1.0, 2.0, 3.0]], dtype=np.float64
    )

    ds = load_unstructured_mesh(connect, coords, data)

    sl = SlicePlot(ds, "z", ("connect1", "test"))
    sl.annotate_mesh_lines()

    sl = SlicePlot(ds, "z", ("connect1", "testAgain"))
    sl.annotate_mesh_lines()

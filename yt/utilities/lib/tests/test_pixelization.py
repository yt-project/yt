import numpy as np

from yt.utilities.lib.pixelization_routines import sample_arbitrary_points_in_1d_buffer


def test_sample_arbitrary_points_in_1d_buffer():
    xbuff = np.array([0.52, 0.54, 0.3], dtype=np.float64)
    ybuff = np.array([0.52, 0.49, 0.3], dtype=np.float64)
    zbuff = np.array([1.0, 0.78, 0.3], dtype=np.float64)
    buff = np.zeros(xbuff.shape)

    ele_0 = np.array([0.5, 0.5, 0.5], dtype=np.float64)
    ele_1 = np.array([0.5, 0.5, 0.5], dtype=np.float64)
    ele_2 = np.array([0.75, 0.98, 1.5], dtype=np.float64)

    dele_0 = np.array([0.1, 0.1, 0.1], dtype=np.float64)
    dele_1 = np.array([0.1, 0.1, 0.1], dtype=np.float64)
    dele_2 = np.array([0.1, 0.1, 0.1], dtype=np.float64)

    data = np.array([0.1, 0.2, 0.3], dtype=np.float64)

    msk = sample_arbitrary_points_in_1d_buffer(
        buff,
        xbuff,
        ybuff,
        zbuff,
        ele_0,
        ele_1,
        ele_2,
        dele_0,
        dele_1,
        dele_2,
        data,
        return_mask=1,
    )
    assert buff[0] == data[1]
    assert buff[1] == data[0]
    assert buff[2] == 0.0
    assert not msk[2]

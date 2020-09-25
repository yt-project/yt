import numpy as np

from yt.utilities.lib.image_samplers import run_interpolation_routines


def test_fast_interpolate():
    np.random.seed(0x4D3D3D3)
    shape = (32, 16, 7)
    N = 10000
    M = 1000
    random_values = np.random.random(shape)
    random_index = np.empty((M, 3), dtype="u8")
    for i in range(3):
        random_index[:, i] = np.random.randint(0, shape[i], size=M, dtype="u8")
    random_offset = np.random.random((M, 3))
    # import time
    # t1 = time.time()
    ngood, nbad = run_interpolation_routines(
        "fast", N, random_values, random_index, random_offset
    )
    # t2 = time.time()
    # print(f"Took {t2-t1} seconds with {nbad} bad results and {ngood} good results")
    # print(random_values.sum(dtype="float64"), random_values.mean(dtype="float64"))
    # print(random_index.sum(dtype="float64"), random_index.mean(dtype="float64"))
    # print(random_offset.sum(dtype="float64"), random_offset.mean(dtype="float64"))
    assert nbad == 0


def test_offset_interpolate():
    pass


def test_trilinear_interpolate():
    pass


def test_eval_gradient():
    pass


if __name__ == "__main__":
    test_fast_interpolate()

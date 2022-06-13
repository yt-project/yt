import numpy as np

from yt.testing import (
    assert_array_equal,
    assert_array_less,
    assert_equal,
    assert_raises,
    fake_random_ds,
)
from yt.utilities.lib.misc_utilities import (
    obtain_position_vector,
    obtain_relative_velocity_vector,
)

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")
_units = ("g/cm**3", "cm/s", "cm/s", "cm/s")

# TODO: error compact/spread bits for incorrect size
# TODO: test msdb for [0,0], [1,1], [2,2] etc.


def test_spread_bits():
    from yt.utilities.lib.geometry_utils import spread_bits

    li = [
        (
            np.uint64(0b111111111111111111111),
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001001),
        )
    ]
    for i, ans in li:
        out = spread_bits(i)
        assert_equal(out, ans)


def test_compact_bits():
    from yt.utilities.lib.geometry_utils import compact_bits

    li = [
        (
            np.uint64(0b111111111111111111111),
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001001),
        )
    ]
    for ans, i in li:
        out = compact_bits(i)
        assert_equal(out, ans)


def test_spread_and_compact_bits():
    from yt.utilities.lib.geometry_utils import compact_bits, spread_bits

    li = [np.uint64(0b111111111111111111111)]
    for ans in li:
        mi = spread_bits(ans)
        out = compact_bits(mi)
        assert_equal(out, ans)


def test_lsz():
    from yt.utilities.lib.geometry_utils import lsz

    li = [
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001001),
            3 * 21,
            3,
            0,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001000),
            3 * 0,
            3,
            0,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001000001),
            3 * 1,
            3,
            0,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001000001001),
            3 * 2,
            3,
            0,
        ),
        (
            np.uint64(0b10010010010010010010010010010010010010010010010010010010010010),
            3 * 0,
            3,
            0,
        ),
        (
            np.uint64(
                0b100100100100100100100100100100100100100100100100100100100100100
            ),
            3 * 0,
            3,
            0,
        ),
        (np.uint64(0b100), 0, 1, 0),
        (np.uint64(0b100), 1, 1, 1),
        (np.uint64(0b100), 3, 1, 2),
        (np.uint64(0b100), 3, 1, 3),
    ]
    for i, ans, stride, start in li:
        out = lsz(i, stride=stride, start=start)
        assert_equal(out, ans)


def test_lsb():
    from yt.utilities.lib.geometry_utils import lsb

    li = [
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001001),
            3 * 0,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001001000),
            3 * 1,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001001000000),
            3 * 2,
        ),
        (
            np.uint64(0b1001001001001001001001001001001001001001001001001001000000000),
            3 * 3,
        ),
        (
            np.uint64(0b10010010010010010010010010010010010010010010010010010010010010),
            3 * 21,
        ),
        (
            np.uint64(
                0b100100100100100100100100100100100100100100100100100100100100100
            ),
            3 * 21,
        ),
    ]
    for i, ans in li:
        out = lsb(i, stride=3)
        assert_equal(out, ans)


def test_bitwise_addition():
    from yt.utilities.lib.geometry_utils import bitwise_addition

    # TODO: Handle negative & periodic boundaries
    lz = [
        (0, 1),
        #          (0,-1),
        (1, 1),
        (1, 2),
        (1, 4),
        (1, -1),
        (2, 1),
        (2, 2),
        (2, -1),
        (2, -2),
        (3, 1),
        (3, 5),
        (3, -1),
    ]
    for i, a in lz:
        i = np.uint64(i)
        a = np.int64(a)
        out = bitwise_addition(i, a, stride=1, start=0)
        assert_equal(out, i + a)


# def test_add_to_morton_coord():
#    from yt.utilities.lib.geometry_utils import add_to_morton_coord


def test_get_morton_indices():
    from yt.utilities.lib.geometry_utils import (
        get_morton_indices,
        get_morton_indices_unravel,
    )

    INDEX_MAX_64 = np.uint64(2097151)
    li = np.arange(6, dtype=np.uint64).reshape((2, 3))
    mi_ans = np.array([10, 229], dtype=np.uint64)
    mi_out = get_morton_indices(li)
    mi_out2 = get_morton_indices_unravel(li[:, 0], li[:, 1], li[:, 2])
    assert_array_equal(mi_out, mi_ans)
    assert_array_equal(mi_out2, mi_ans)
    li[0, :] = INDEX_MAX_64 * np.ones(3, dtype=np.uint64)
    assert_raises(ValueError, get_morton_indices, li)
    assert_raises(ValueError, get_morton_indices_unravel, li[:, 0], li[:, 1], li[:, 2])


def test_get_morton_points():
    from yt.utilities.lib.geometry_utils import get_morton_points

    mi = np.array([10, 229], dtype=np.uint64)
    li_ans = np.arange(6, dtype=np.uint64).reshape((2, 3))
    li_out = get_morton_points(mi)
    assert_array_equal(li_out, li_ans)


def test_compare_morton():
    # TODO: Add error messages to assertions
    from yt.utilities.lib.geometry_utils import compare_morton

    # Diagonal
    p = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    assert_equal(compare_morton(p, q), 1)
    assert_equal(compare_morton(q, p), 0)
    assert_equal(compare_morton(p, p), 0)
    # 1-1 vs 0-1
    p = np.array([1.0, 1.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    assert_equal(compare_morton(p, q), 1)
    assert_equal(compare_morton(q, p), 0)
    assert_equal(compare_morton(p, p), 0)
    # x advance, y decrease
    p = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    assert_equal(compare_morton(p, q), 1)
    assert_equal(compare_morton(q, p), 0)
    assert_equal(compare_morton(p, p), 0)
    # x&y advance, z decrease
    p = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    q = np.array([1.0, 1.0, 0.0], dtype=np.float64)
    assert_equal(compare_morton(p, q), 1)
    assert_equal(compare_morton(q, p), 0)
    assert_equal(compare_morton(p, p), 0)


def test_get_morton_neighbors_coarse():
    from yt.utilities.lib.geometry_utils import get_morton_neighbors_coarse

    imax = 5
    ngz = 1
    tests = {
        (7, 1): np.array(
            [
                35,
                49,
                56,
                48,
                33,
                40,
                32,
                42,
                34,
                3,
                17,
                24,
                16,
                1,
                8,
                0,
                10,
                2,
                21,
                28,
                20,
                5,
                12,
                4,
                14,
                6,
            ],
            dtype="uint64",
        ),
        (7, 0): np.array(
            [
                35,
                49,
                56,
                48,
                33,
                40,
                32,
                42,
                34,
                3,
                17,
                24,
                16,
                1,
                8,
                0,
                10,
                2,
                21,
                28,
                20,
                5,
                12,
                4,
                14,
                6,
            ],
            dtype="uint64",
        ),
        (0, 1): np.array(
            [
                4,
                6,
                7,
                70,
                132,
                133,
                196,
                5,
                68,
                256,
                258,
                259,
                322,
                384,
                385,
                448,
                257,
                320,
                2,
                3,
                66,
                128,
                129,
                192,
                1,
                64,
            ],
            dtype="uint64",
        ),
        (0, 0): np.array([4, 6, 7, 5, 2, 3, 1], dtype="uint64"),
        (448, 1): np.array(
            [
                192,
                64,
                0,
                9,
                82,
                18,
                27,
                128,
                137,
                228,
                100,
                36,
                45,
                118,
                54,
                63,
                164,
                173,
                320,
                256,
                265,
                338,
                274,
                283,
                384,
                393,
            ],
            dtype="uint64",
        ),
        (448, 0): np.array([228, 118, 63, 173, 338, 283, 393], dtype="uint64"),
    }
    for (mi1, periodic), ans in tests.items():
        n1 = get_morton_neighbors_coarse(mi1, imax, periodic, ngz)
        assert_equal(np.sort(n1), np.sort(ans))


def test_get_morton_neighbors_refined():
    from yt.utilities.lib.geometry_utils import get_morton_neighbors_refined

    imax1 = 5
    imax2 = 5
    ngz = 1
    tests = {
        (7, 7, 1): (
            np.array(
                [
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                ],
                dtype="uint64",
            ),
            np.array(
                [
                    35,
                    49,
                    56,
                    48,
                    33,
                    40,
                    32,
                    42,
                    34,
                    3,
                    17,
                    24,
                    16,
                    1,
                    8,
                    0,
                    10,
                    2,
                    21,
                    28,
                    20,
                    5,
                    12,
                    4,
                    14,
                    6,
                ],
                dtype="uint64",
            ),
        ),
        (7, 7, 0): (
            np.array(
                [
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                    7,
                ],
                dtype="uint64",
            ),
            np.array(
                [
                    35,
                    49,
                    56,
                    48,
                    33,
                    40,
                    32,
                    42,
                    34,
                    3,
                    17,
                    24,
                    16,
                    1,
                    8,
                    0,
                    10,
                    2,
                    21,
                    28,
                    20,
                    5,
                    12,
                    4,
                    14,
                    6,
                ],
                dtype="uint64",
            ),
        ),
        (0, 0, 1): (
            np.array(
                [
                    0,
                    0,
                    0,
                    64,
                    128,
                    128,
                    192,
                    0,
                    64,
                    256,
                    256,
                    256,
                    320,
                    384,
                    384,
                    448,
                    256,
                    320,
                    0,
                    0,
                    64,
                    128,
                    128,
                    192,
                    0,
                    64,
                ],
                dtype="uint64",
            ),
            np.array(
                [
                    4,
                    6,
                    7,
                    70,
                    132,
                    133,
                    196,
                    5,
                    68,
                    256,
                    258,
                    259,
                    322,
                    384,
                    385,
                    448,
                    257,
                    320,
                    2,
                    3,
                    66,
                    128,
                    129,
                    192,
                    1,
                    64,
                ],
                dtype="uint64",
            ),
        ),
        (0, 0, 0): (
            np.array([0, 0, 0, 0, 0, 0, 0], dtype="uint64"),
            np.array([4, 6, 7, 5, 2, 3, 1], dtype="uint64"),
        ),
        (448, 448, 1): (
            np.array(
                [
                    192,
                    64,
                    0,
                    64,
                    192,
                    128,
                    192,
                    128,
                    192,
                    448,
                    320,
                    256,
                    320,
                    448,
                    384,
                    448,
                    384,
                    448,
                    320,
                    256,
                    320,
                    448,
                    384,
                    448,
                    384,
                    448,
                ],
                dtype="uint64",
            ),
            np.array(
                [
                    192,
                    64,
                    0,
                    9,
                    82,
                    18,
                    27,
                    128,
                    137,
                    228,
                    100,
                    36,
                    45,
                    118,
                    54,
                    63,
                    164,
                    173,
                    320,
                    256,
                    265,
                    338,
                    274,
                    283,
                    384,
                    393,
                ],
                dtype="uint64",
            ),
        ),
        (448, 448, 0): (
            np.array([448, 448, 448, 448, 448, 448, 448], dtype="uint64"),
            np.array([228, 118, 63, 173, 338, 283, 393], dtype="uint64"),
        ),
    }
    for (mi1, mi2, periodic), (ans1, ans2) in tests.items():
        n1, n2 = get_morton_neighbors_refined(mi1, mi2, imax1, imax2, periodic, ngz)
        assert_equal(np.sort(n1), np.sort(ans1))
        assert_equal(np.sort(n2), np.sort(ans2))


def test_morton_neighbor():
    from yt.utilities.lib.geometry_utils import get_morton_indices, morton_neighbor

    order = 20
    imax = np.uint64(1 << order)
    p = np.array(
        [
            [imax / 2, imax / 2, imax / 2],
            [imax / 2, imax / 2, 0],
            [imax / 2, imax / 2, imax],
        ],
        dtype=np.uint64,
    )
    p_ans = np.array(
        [
            [imax / 2, imax / 2, imax / 2 + 1],
            [imax / 2, imax / 2, imax / 2 - 1],
            [imax / 2, imax / 2, imax - 1],
            [imax / 2, imax / 2, 1],
            [imax / 2, imax / 2 + 1, imax / 2 + 1],
            [imax / 2 - 1, imax / 2 - 1, imax / 2],
            [imax / 2 - 1, imax / 2, imax / 2 + 1],
            [imax / 2, imax / 2 - 1, imax - 1],
            [imax / 2, imax / 2 + 1, 1],
        ],
        dtype=np.uint64,
    )
    mi_ans = get_morton_indices(p_ans)
    assert_equal(morton_neighbor(p[0, :], [2], [+1], imax), mi_ans[0])
    assert_equal(morton_neighbor(p[0, :], [2], [-1], imax), mi_ans[1])
    assert_equal(morton_neighbor(p[1, :], [2], [-1], imax, periodic=False), -1)
    assert_equal(morton_neighbor(p[2, :], [2], [+1], imax, periodic=False), -1)
    assert_equal(morton_neighbor(p[1, :], [2], [-1], imax, periodic=True), mi_ans[2])
    assert_equal(morton_neighbor(p[2, :], [2], [+1], imax, periodic=True), mi_ans[3])
    assert_equal(morton_neighbor(p[0, :], [1, 2], [+1, +1], imax), mi_ans[4])
    assert_equal(morton_neighbor(p[0, :], [0, 1], [-1, -1], imax), mi_ans[5])
    assert_equal(morton_neighbor(p[0, :], [0, 2], [-1, +1], imax), mi_ans[6])
    assert_equal(morton_neighbor(p[1, :], [1, 2], [-1, -1], imax, periodic=False), -1)
    assert_equal(morton_neighbor(p[2, :], [1, 2], [+1, +1], imax, periodic=False), -1)
    assert_equal(
        morton_neighbor(p[1, :], [1, 2], [-1, -1], imax, periodic=True), mi_ans[7]
    )
    assert_equal(
        morton_neighbor(p[2, :], [1, 2], [+1, +1], imax, periodic=True), mi_ans[8]
    )


def test_get_morton_neighbors():
    from yt.utilities.lib.geometry_utils import get_morton_indices, get_morton_neighbors

    order = 20
    imax = 1 << order
    p = np.array(
        [
            [imax / 2, imax / 2, imax / 2],
            [imax / 2, imax / 2, 0],
            [imax / 2, imax / 2, imax],
        ],
        dtype=np.uint64,
    )
    pn_non = [
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, imax / 2],
                [imax / 2 + 1, imax / 2 + 1, imax / 2],
                [imax / 2 + 1, imax / 2 + 1, imax / 2 + 1],
                [imax / 2 + 1, imax / 2 + 1, imax / 2 - 1],
                [imax / 2 + 1, imax / 2 - 1, imax / 2],
                [imax / 2 + 1, imax / 2 - 1, imax / 2 + 1],
                [imax / 2 + 1, imax / 2 - 1, imax / 2 - 1],
                [imax / 2 + 1, imax / 2, imax / 2 + 1],
                [imax / 2 + 1, imax / 2, imax / 2 - 1],
                [imax / 2 - 1, imax / 2, imax / 2],
                [imax / 2 - 1, imax / 2 + 1, imax / 2],
                [imax / 2 - 1, imax / 2 + 1, imax / 2 + 1],
                [imax / 2 - 1, imax / 2 + 1, imax / 2 - 1],
                [imax / 2 - 1, imax / 2 - 1, imax / 2],
                [imax / 2 - 1, imax / 2 - 1, imax / 2 + 1],
                [imax / 2 - 1, imax / 2 - 1, imax / 2 - 1],
                [imax / 2 - 1, imax / 2, imax / 2 + 1],
                [imax / 2 - 1, imax / 2, imax / 2 - 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, imax / 2],
                [imax / 2, imax / 2 + 1, imax / 2 + 1],
                [imax / 2, imax / 2 + 1, imax / 2 - 1],
                [imax / 2, imax / 2 - 1, imax / 2],
                [imax / 2, imax / 2 - 1, imax / 2 + 1],
                [imax / 2, imax / 2 - 1, imax / 2 - 1],
                # x +/- 1
                [imax / 2, imax / 2, imax / 2 + 1],
                [imax / 2, imax / 2, imax / 2 - 1],
            ],
            dtype=np.uint64,
        ),
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, 0],
                [imax / 2 + 1, imax / 2 + 1, 0],
                [imax / 2 + 1, imax / 2 + 1, 1],
                [imax / 2 + 1, imax / 2 - 1, 0],
                [imax / 2 + 1, imax / 2 - 1, 1],
                [imax / 2 + 1, imax / 2, 1],
                [imax / 2 - 1, imax / 2, 0],
                [imax / 2 - 1, imax / 2 + 1, 0],
                [imax / 2 - 1, imax / 2 + 1, 1],
                [imax / 2 - 1, imax / 2 - 1, 0],
                [imax / 2 - 1, imax / 2 - 1, 1],
                [imax / 2 - 1, imax / 2, 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, 0],
                [imax / 2, imax / 2 + 1, 1],
                [imax / 2, imax / 2 - 1, 0],
                [imax / 2, imax / 2 - 1, 1],
                # z +/- 1
                [imax / 2, imax / 2, 0 + 1],
            ],
            dtype=np.uint64,
        ),
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, imax],
                [imax / 2 + 1, imax / 2 + 1, imax],
                [imax / 2 + 1, imax / 2 + 1, imax - 1],
                [imax / 2 + 1, imax / 2 - 1, imax],
                [imax / 2 + 1, imax / 2 - 1, imax - 1],
                [imax / 2 + 1, imax / 2, imax - 1],
                [imax / 2 - 1, imax / 2, imax],
                [imax / 2 - 1, imax / 2 + 1, imax],
                [imax / 2 - 1, imax / 2 + 1, imax - 1],
                [imax / 2 - 1, imax / 2 - 1, imax],
                [imax / 2 - 1, imax / 2 - 1, imax - 1],
                [imax / 2 - 1, imax / 2, imax - 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, imax],
                [imax / 2, imax / 2 + 1, imax - 1],
                [imax / 2, imax / 2 - 1, imax],
                [imax / 2, imax / 2 - 1, imax - 1],
                # z +/- 1
                [imax / 2, imax / 2, imax - 1],
            ],
            dtype=np.uint64,
        ),
    ]
    pn_per = [
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, imax / 2],
                [imax / 2 + 1, imax / 2 + 1, imax / 2],
                [imax / 2 + 1, imax / 2 + 1, imax / 2 + 1],
                [imax / 2 + 1, imax / 2 + 1, imax / 2 - 1],
                [imax / 2 + 1, imax / 2 - 1, imax / 2],
                [imax / 2 + 1, imax / 2 - 1, imax / 2 + 1],
                [imax / 2 + 1, imax / 2 - 1, imax / 2 - 1],
                [imax / 2 + 1, imax / 2, imax / 2 + 1],
                [imax / 2 + 1, imax / 2, imax / 2 - 1],
                [imax / 2 - 1, imax / 2, imax / 2],
                [imax / 2 - 1, imax / 2 + 1, imax / 2],
                [imax / 2 - 1, imax / 2 + 1, imax / 2 + 1],
                [imax / 2 - 1, imax / 2 + 1, imax / 2 - 1],
                [imax / 2 - 1, imax / 2 - 1, imax / 2],
                [imax / 2 - 1, imax / 2 - 1, imax / 2 + 1],
                [imax / 2 - 1, imax / 2 - 1, imax / 2 - 1],
                [imax / 2 - 1, imax / 2, imax / 2 + 1],
                [imax / 2 - 1, imax / 2, imax / 2 - 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, imax / 2],
                [imax / 2, imax / 2 + 1, imax / 2 + 1],
                [imax / 2, imax / 2 + 1, imax / 2 - 1],
                [imax / 2, imax / 2 - 1, imax / 2],
                [imax / 2, imax / 2 - 1, imax / 2 + 1],
                [imax / 2, imax / 2 - 1, imax / 2 - 1],
                # z +/- 1
                [imax / 2, imax / 2, imax / 2 + 1],
                [imax / 2, imax / 2, imax / 2 - 1],
            ],
            dtype=np.uint64,
        ),
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, 0],
                [imax / 2 + 1, imax / 2 + 1, 0],
                [imax / 2 + 1, imax / 2 + 1, 1],
                [imax / 2 + 1, imax / 2 + 1, imax - 1],
                [imax / 2 + 1, imax / 2 - 1, 0],
                [imax / 2 + 1, imax / 2 - 1, 1],
                [imax / 2 + 1, imax / 2 - 1, imax - 1],
                [imax / 2 + 1, imax / 2, 1],
                [imax / 2 + 1, imax / 2, imax - 1],
                [imax / 2 - 1, imax / 2, 0],
                [imax / 2 - 1, imax / 2 + 1, 0],
                [imax / 2 - 1, imax / 2 + 1, 1],
                [imax / 2 - 1, imax / 2 + 1, imax - 1],
                [imax / 2 - 1, imax / 2 - 1, 0],
                [imax / 2 - 1, imax / 2 - 1, 1],
                [imax / 2 - 1, imax / 2 - 1, imax - 1],
                [imax / 2 - 1, imax / 2, 1],
                [imax / 2 - 1, imax / 2, imax - 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, 0],
                [imax / 2, imax / 2 + 1, 1],
                [imax / 2, imax / 2 + 1, imax - 1],
                [imax / 2, imax / 2 - 1, 0],
                [imax / 2, imax / 2 - 1, 1],
                [imax / 2, imax / 2 - 1, imax - 1],
                # z +/- 1
                [imax / 2, imax / 2, 0 + 1],
                [imax / 2, imax / 2, imax - 1],
            ],
            dtype=np.uint64,
        ),
        np.array(
            [
                # x +/- 1
                [imax / 2 + 1, imax / 2, imax],
                [imax / 2 + 1, imax / 2 + 1, imax],
                [imax / 2 + 1, imax / 2 + 1, 1],
                [imax / 2 + 1, imax / 2 + 1, imax - 1],
                [imax / 2 + 1, imax / 2 - 1, imax],
                [imax / 2 + 1, imax / 2 - 1, 1],
                [imax / 2 + 1, imax / 2 - 1, imax - 1],
                [imax / 2 + 1, imax / 2, 1],
                [imax / 2 + 1, imax / 2, imax - 1],
                [imax / 2 - 1, imax / 2, imax],
                [imax / 2 - 1, imax / 2 + 1, imax],
                [imax / 2 - 1, imax / 2 + 1, 1],
                [imax / 2 - 1, imax / 2 + 1, imax - 1],
                [imax / 2 - 1, imax / 2 - 1, imax],
                [imax / 2 - 1, imax / 2 - 1, 1],
                [imax / 2 - 1, imax / 2 - 1, imax - 1],
                [imax / 2 - 1, imax / 2, 1],
                [imax / 2 - 1, imax / 2, imax - 1],
                # y +/- 1
                [imax / 2, imax / 2 + 1, imax],
                [imax / 2, imax / 2 + 1, 1],
                [imax / 2, imax / 2 + 1, imax - 1],
                [imax / 2, imax / 2 - 1, imax],
                [imax / 2, imax / 2 - 1, 1],
                [imax / 2, imax / 2 - 1, imax - 1],
                # z +/- 1
                [imax / 2, imax / 2, 1],
                [imax / 2, imax / 2, imax - 1],
            ],
            dtype=np.uint64,
        ),
    ]
    mi = get_morton_indices(p)
    N = mi.shape[0]
    # Non-periodic
    for i in range(N):
        out = get_morton_neighbors(
            np.array([mi[i]], dtype=np.uint64), order=order, periodic=False
        )
        ans = get_morton_indices(np.vstack([p[i, :], pn_non[i]]))
        assert_array_equal(np.unique(out), np.unique(ans), err_msg=f"Non-periodic: {i}")
    # Periodic
    for i in range(N):
        out = get_morton_neighbors(
            np.array([mi[i]], dtype=np.uint64), order=order, periodic=True
        )
        ans = get_morton_indices(np.vstack([p[i, :], pn_per[i]]))
        assert_array_equal(np.unique(out), np.unique(ans), err_msg=f"Periodic: {i}")


def test_dist():
    from yt.utilities.lib.geometry_utils import dist

    p = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    q = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    assert_equal(dist(p, q), 0.0)
    p = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    assert_equal(dist(p, q), 1.0)
    p = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 1.0, 0.0], dtype=np.float64)
    assert_equal(dist(p, q), np.sqrt(2.0))
    p = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    q = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    assert_equal(dist(p, q), np.sqrt(3.0))


def test_knn_direct(seed=1):
    from yt.utilities.lib.geometry_utils import knn_direct

    np.random.seed(seed)
    k = 64
    N = 1e5
    idx = np.arange(N, dtype=np.uint64)
    rad = np.arange(N, dtype=np.float64)
    pos = np.vstack(3 * [rad**2 / 3.0]).T
    sort_shf = np.arange(N, dtype=np.uint64)
    for _ in range(20):
        np.random.shuffle(sort_shf)
        sort_ans = np.argsort(sort_shf)[:k]
        sort_out = knn_direct(pos[sort_shf, :], k, sort_ans[0], idx)
        assert_array_equal(sort_out, sort_ans)


# TODO: test of quadtree (.pxd)


def test_obtain_position_vector():
    ds = fake_random_ds(
        64, nprocs=8, fields=_fields, units=_units, negative=[False, True, True, True]
    )

    dd = ds.sphere((0.5, 0.5, 0.5), 0.2)

    coords = obtain_position_vector(dd)

    r = np.sqrt(np.sum(coords * coords, axis=0))

    assert_array_less(r.max(), 0.2)

    assert_array_less(0.0, r.min())


def test_obtain_relative_velocity_vector():
    ds = fake_random_ds(
        64, nprocs=8, fields=_fields, units=_units, negative=[False, True, True, True]
    )

    dd = ds.all_data()

    vels = obtain_relative_velocity_vector(dd)

    assert_array_equal(vels[0, :], dd[("gas", "velocity_x")])
    assert_array_equal(vels[1, :], dd[("gas", "velocity_y")])
    assert_array_equal(vels[2, :], dd[("gas", "velocity_z")])

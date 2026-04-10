import tempfile
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_raises

import yt.utilities.linear_interpolators as lin


def test_unilinear_interpolator_exact_linear_examples():
    interpolator = lin.UnilinearFieldInterpolator(
        np.array([0.0, 10.0]), (0.0, 1.0), "x", truncate=True
    )

    assert_allclose(interpolator({"x": np.array([0.25, 0.5, 0.75])}), [2.5, 5.0, 7.5])
    assert_allclose(interpolator({"x": np.array([-0.5, 1.5])}), [-5.0, 15.0])

    strict_interpolator = lin.UnilinearFieldInterpolator(
        np.array([0.0, 10.0]), (0.0, 1.0), "x", truncate=False
    )
    with assert_raises(ValueError):
        strict_interpolator({"x": np.array([-0.1])})


def test_unilinear_constructor_rejects_mismatched_bins():
    with assert_raises(ValueError):
        lin.UnilinearFieldInterpolator(
            np.array([0.0, 10.0]), np.array([0.0]), "x", truncate=True
        )


def test_bilinear_interpolator_exact_linear_examples():
    interpolator = lin.BilinearFieldInterpolator(
        np.array([[0.0, 2.0], [1.0, 3.0]]),
        (0.0, 1.0, 0.0, 1.0),
        ["x", "y"],
        truncate=True,
    )

    assert_allclose(
        interpolator({"x": np.array([0.25, 0.5]), "y": np.array([0.5, 0.25])}),
        [1.25, 1.0],
    )
    assert_allclose(interpolator({"x": np.array([-0.5]), "y": np.array([0.5])}), [0.5])

    strict_interpolator = lin.BilinearFieldInterpolator(
        np.array([[0.0, 2.0], [1.0, 3.0]]),
        (0.0, 1.0, 0.0, 1.0),
        ["x", "y"],
        truncate=False,
    )
    with assert_raises(ValueError):
        strict_interpolator({"x": np.array([-0.1]), "y": np.array([0.5])})


def test_bilinear_constructor_rejects_invalid_boundary_specifications():
    table = np.array([[0.0, 2.0], [1.0, 3.0]])

    with assert_raises(ValueError):
        lin.BilinearFieldInterpolator(
            table,
            (np.array([0.0]), np.array([0.0, 1.0])),
            ["x", "y"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.BilinearFieldInterpolator(
            table,
            (np.array([0.0, 1.0]), np.array([0.0])),
            ["x", "y"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.BilinearFieldInterpolator(table, (0.0, 1.0, 0.0), ["x", "y"], True)


def test_trilinear_interpolator_exact_linear_examples():
    table = np.zeros((2, 2, 2))
    for i, x_val in enumerate([0.0, 1.0]):
        for j, y_val in enumerate([0.0, 1.0]):
            for k, z_val in enumerate([0.0, 1.0]):
                table[i, j, k] = x_val + 2.0 * y_val + 3.0 * z_val

    interpolator = lin.TrilinearFieldInterpolator(
        table,
        (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
        ["x", "y", "z"],
        truncate=True,
    )

    assert_allclose(
        interpolator(
            {
                "x": np.array([0.25]),
                "y": np.array([0.5]),
                "z": np.array([0.75]),
            }
        ),
        [3.5],
    )
    assert_allclose(
        interpolator(
            {
                "x": np.array([0.25, -0.5]),
                "y": np.array([0.5, 0.5]),
                "z": np.array([0.75, 0.75]),
            }
        ),
        [3.5, 2.75],
    )

    strict_interpolator = lin.TrilinearFieldInterpolator(
        table,
        (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
        ["x", "y", "z"],
        truncate=False,
    )
    with assert_raises(ValueError):
        strict_interpolator(
            {"x": np.array([-0.1]), "y": np.array([0.5]), "z": np.array([0.75])}
        )


def test_trilinear_constructor_and_strict_mode_cover_y_and_z_errors():
    table = np.zeros((2, 2, 2))

    with assert_raises(ValueError):
        lin.TrilinearFieldInterpolator(
            table,
            (np.array([0.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0])),
            ["x", "y", "z"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.TrilinearFieldInterpolator(
            table,
            (np.array([0.0, 1.0]), np.array([0.0]), np.array([0.0, 1.0])),
            ["x", "y", "z"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.TrilinearFieldInterpolator(
            table,
            (np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([0.0])),
            ["x", "y", "z"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.TrilinearFieldInterpolator(
            table,
            (np.array([0.0, 1.0]), np.array([0.0, 1.0])),
            ["x", "y", "z"],
            truncate=True,
        )

    strict_interpolator = lin.TrilinearFieldInterpolator(
        table,
        (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
        ["x", "y", "z"],
        truncate=False,
    )
    with assert_raises(ValueError):
        strict_interpolator(
            {"x": np.array([0.5]), "y": np.array([1.1]), "z": np.array([0.5])}
        )
    with assert_raises(ValueError):
        strict_interpolator(
            {"x": np.array([0.5]), "y": np.array([0.5]), "z": np.array([1.1])}
        )


def test_quadrilinear_interpolator_exact_linear_examples():
    table = np.zeros((2, 2, 2, 2))
    for i, x_val in enumerate([0.0, 1.0]):
        for j, y_val in enumerate([0.0, 1.0]):
            for k, z_val in enumerate([0.0, 1.0]):
                for l, w_val in enumerate([0.0, 1.0]):
                    table[i, j, k, l] = x_val + 2.0 * y_val + 3.0 * z_val + 4.0 * w_val

    interpolator = lin.QuadrilinearFieldInterpolator(
        table,
        (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
        ["x", "y", "z", "w"],
        truncate=True,
    )
    assert_allclose(
        interpolator(
            {
                "x": np.array([0.25, -0.5]),
                "y": np.array([0.5, 0.5]),
                "z": np.array([0.75, 0.75]),
                "w": np.array([0.125, 0.125]),
            }
        ),
        [4.0, 3.25],
    )

    array_boundary_interpolator = lin.QuadrilinearFieldInterpolator(
        table,
        (
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
        ),
        ["x", "y", "z", "w"],
        truncate=True,
    )
    assert_allclose(
        array_boundary_interpolator(
            {
                "x": np.array([0.25]),
                "y": np.array([0.5]),
                "z": np.array([0.75]),
                "w": np.array([0.125]),
            }
        ),
        [4.0],
    )

    strict_interpolator = lin.QuadrilinearFieldInterpolator(
        table,
        (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
        ["x", "y", "z", "w"],
        truncate=False,
    )
    with assert_raises(ValueError):
        strict_interpolator(
            {
                "x": np.array([0.5]),
                "y": np.array([0.5]),
                "z": np.array([0.5]),
                "w": np.array([1.1]),
            }
        )


def test_quadrilinear_constructor_rejects_invalid_boundary_specifications():
    table = np.zeros((2, 2, 2, 2))

    with assert_raises(ValueError):
        lin.QuadrilinearFieldInterpolator(
            table,
            (
                np.array([0.0]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
            ),
            ["x", "y", "z", "w"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.QuadrilinearFieldInterpolator(
            table,
            (
                np.array([0.0, 1.0]),
                np.array([0.0]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
            ),
            ["x", "y", "z", "w"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.QuadrilinearFieldInterpolator(
            table,
            (
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([0.0]),
                np.array([0.0, 1.0]),
            ),
            ["x", "y", "z", "w"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.QuadrilinearFieldInterpolator(
            table,
            (
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([0.0]),
            ),
            ["x", "y", "z", "w"],
            truncate=True,
        )
    with assert_raises(ValueError):
        lin.QuadrilinearFieldInterpolator(
            table,
            (
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
            ),
            ["x", "y", "z", "w"],
            truncate=True,
        )


class _FakeSphereDataset:
    def __init__(self, unit_scale):
        self.unit_scale = unit_scale

    def __getitem__(self, item):
        assert item == "code_length"
        return self.unit_scale

    def sphere(self, center, radius):
        return center, radius


def test_get_centers_skips_comments_and_scales_radius_by_unit():
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = Path(tmpdir) / "centers.dat"
        filename.write_text(
            "# x y z r\n"
            "0.1 0.2 0.3 4.0\n"
            "0.4 0.5 0.6 2.0\n",
            encoding="ascii",
        )

        centers = list(
            lin.get_centers(
                _FakeSphereDataset(2.0),
                str(filename),
                center_cols=(0, 1, 2),
                radius_col=3,
                unit="code_length",
            )
        )

    assert centers == [
        ([0.1, 0.2, 0.3], 2.0),
        ([0.4, 0.5, 0.6], 1.0),
    ]

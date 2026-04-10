import numpy as np
from numpy.testing import assert_allclose

import yt.utilities.initial_conditions as ic


class _ArrayWrapper:
    def __init__(self, values):
        self.values = np.array(values, dtype="float64")

    def to_ndarray(self):
        return self.values

    def ndarray_view(self):
        return self.values

    def __sub__(self, other):
        return self.values - other

    def __rsub__(self, other):
        return other - self.values

    def __getitem__(self, item):
        return self.values[item]

    def __setitem__(self, item, value):
        self.values[item] = value


class _FakeGrid(dict):
    def __init__(self, x, y, z, **fields):
        super().__init__()
        self.ActiveDimensions = np.array(x).shape
        self["x"] = _ArrayWrapper(x)
        self["y"] = _ArrayWrapper(y)
        self["z"] = _ArrayWrapper(z)
        for field, values in fields.items():
            self[field] = np.array(values, dtype="float64")


class _FakeIndex:
    def __init__(self, grids):
        self.grids = grids


class _FakeDS:
    def __init__(self, grids):
        self.index = _FakeIndex(grids)


class _FixedRNG:
    def __init__(self, values):
        self.values = np.array(values, dtype="float64")

    def random(self, shape):
        assert shape == self.values.shape
        return self.values


class _SequenceRNG:
    def __init__(self, values):
        self.values = [np.array(value, dtype="float64") for value in values]

    def random(self, shape):
        value = self.values.pop(0)
        assert shape == value.shape
        return value


def test_top_hat_sphere_apply_and_subselection():
    x = np.array([[[0.25]], [[0.50]], [[0.75]], [[0.90]]])
    y = np.array([[[0.50]], [[0.50]], [[0.50]], [[0.50]]])
    z = np.array([[[0.50]], [[0.50]], [[0.50]], [[0.50]]])

    grid = _FakeGrid(x, y, z, density=np.zeros_like(x))
    top_hat = ic.TopHatSphere(0.30, np.array([0.50, 0.50, 0.50]), {"density": 7.0})

    top_hat(grid, sub_select=np.array([[[True]], [[True]], [[False]], [[True]]]))
    assert_allclose(grid["density"][:, 0, 0], [7.0, 7.0, 0.0, 0.0])

    grid_a = _FakeGrid(x, y, z, density=np.zeros_like(x))
    grid_b = _FakeGrid(x, y, z, density=np.zeros_like(x))
    top_hat.apply(_FakeDS([grid_a, grid_b]))
    assert_allclose(grid_a["density"][:, 0, 0], [7.0, 7.0, 7.0, 0.0])
    assert_allclose(grid_b["density"][:, 0, 0], [7.0, 7.0, 7.0, 0.0])


def test_top_hat_sphere_updates_multiple_fields_with_exact_outputs():
    x = np.array([[[0.25]], [[0.50]], [[0.75]], [[0.90]]])
    y = np.array([[[0.50]], [[0.50]], [[0.50]], [[0.50]]])
    z = np.array([[[0.50]], [[0.50]], [[0.50]], [[0.50]]])

    grid = _FakeGrid(
        x,
        y,
        z,
        density=np.zeros_like(x),
        temperature=np.full_like(x, -10.0),
    )
    top_hat = ic.TopHatSphere(
        0.30,
        np.array([0.50, 0.50, 0.50]),
        {"density": 7.0, "temperature": 300.0},
    )

    top_hat(grid)

    assert_allclose(grid["density"][:, 0, 0], [7.0, 7.0, 7.0, 0.0])
    assert_allclose(grid["temperature"][:, 0, 0], [300.0, 300.0, 300.0, -10.0])


def test_cored_sphere_profile_hits_inner_and_outer_values():
    x = np.array([[[0.00]], [[0.50]], [[0.75]], [[1.00]], [[1.25]]])
    zeros = np.zeros_like(x)
    grid = _FakeGrid(x, zeros, zeros, density=np.full_like(x, -1.0))

    cored = ic.CoredSphere(
        core_radius=0.5,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": (10.0, 2.0)},
    )
    cored(grid)

    assert_allclose(
        grid["density"][:, 0, 0],
        [2.0, 2.0, 7.163977794943222, 10.0, -1.0],
    )


def test_beta_model_sphere_matches_expected_profile():
    x = np.array([[[0.00]], [[0.50]], [[0.75]], [[1.00]], [[1.25]]])
    zeros = np.zeros_like(x)
    grid = _FakeGrid(x, zeros, zeros, density=np.zeros_like(x))

    beta_model = ic.BetaModelSphere(
        beta=2.0 / 3.0,
        core_radius=0.5,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": 10.0},
    )
    beta_model(grid)

    assert_allclose(
        grid["density"][:, 0, 0],
        [10.0, 5.0, 3.076923076923077, 2.0, 0.0],
    )


def test_nfw_model_sphere_matches_expected_profile():
    x = np.array([[[0.50]], [[0.75]], [[1.00]], [[1.25]]])
    zeros = np.zeros_like(x)
    grid = _FakeGrid(x, zeros, zeros, density=np.zeros_like(x))

    nfw = ic.NFWModelSphere(
        scale_radius=1.0,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": 8.0},
    )
    nfw(grid)

    assert_allclose(
        grid["density"][:, 0, 0],
        [7.111111111111111, 3.482993197278912, 2.0, 0.0],
    )


def test_random_fluctuation_uses_rng_scale_factor():
    grid = _FakeGrid(
        np.zeros((2, 2, 1)),
        np.zeros((2, 2, 1)),
        np.zeros((2, 2, 1)),
        density=np.array([[[1.0], [2.0]], [[3.0], [4.0]]]),
    )
    operator = ic.RandomFluctuation({"density": 0.2})

    original_default_rng = ic.np.random.default_rng
    try:
        ic.np.random.default_rng = lambda: _FixedRNG(
            np.array([[[0.0], [0.5]], [[1.0], [0.25]]])
        )
        operator(grid)
    finally:
        ic.np.random.default_rng = original_default_rng

    expected = np.array([[[0.9], [2.0]], [[3.3], [3.8]]])
    assert_allclose(grid["density"], expected)


def test_random_fluctuation_updates_multiple_fields_with_exact_scalings():
    grid = _FakeGrid(
        np.zeros((2, 1, 1)),
        np.zeros((2, 1, 1)),
        np.zeros((2, 1, 1)),
        density=np.array([[[10.0]], [[20.0]]]),
        temperature=np.array([[[100.0]], [[200.0]]]),
    )
    operator = ic.RandomFluctuation({"density": 0.2, "temperature": 0.4})

    original_default_rng = ic.np.random.default_rng
    try:
        ic.np.random.default_rng = lambda: _SequenceRNG(
            [
                np.array([[[0.0]], [[1.0]]]),
                np.array([[[0.25]], [[0.75]]]),
            ]
        )
        operator(grid)
    finally:
        ic.np.random.default_rng = original_default_rng

    assert_allclose(grid["density"][:, 0, 0], [9.0, 22.0])
    assert_allclose(grid["temperature"][:, 0, 0], [90.0, 220.0])


def test_profile_operators_and_random_fluctuation_respect_subselection():
    mask = np.array([[[True]], [[False]]])
    zeros = np.zeros((2, 1, 1))
    x_full = np.array([[[0.0]], [[1.0]]])

    cored_grid = _FakeGrid(x_full, zeros, zeros, density=np.full_like(x_full, -1.0))
    ic.CoredSphere(
        core_radius=0.5,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": (10.0, 2.0)},
    )(cored_grid, sub_select=mask)
    assert_allclose(cored_grid["density"][:, 0, 0], [2.0, -1.0])

    beta_grid = _FakeGrid(x_full, zeros, zeros, density=np.zeros_like(x_full))
    ic.BetaModelSphere(
        beta=2.0 / 3.0,
        core_radius=0.5,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": 10.0},
    )(beta_grid, sub_select=mask)
    assert_allclose(beta_grid["density"][:, 0, 0], [10.0, 0.0])

    nfw_grid = _FakeGrid(
        np.array([[[0.5]], [[1.0]]]),
        zeros,
        zeros,
        density=np.zeros_like(x_full),
    )
    ic.NFWModelSphere(
        scale_radius=1.0,
        radius=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        fields={"density": 8.0},
    )(nfw_grid, sub_select=mask)
    assert_allclose(nfw_grid["density"][:, 0, 0], [8.0 / 1.125, 0.0])

    random_grid = _FakeGrid(
        np.zeros((2, 1, 1)),
        np.zeros((2, 1, 1)),
        np.zeros((2, 1, 1)),
        density=np.array([[[1.0]], [[4.0]]]),
    )
    operator = ic.RandomFluctuation({"density": 0.2})

    original_default_rng = ic.np.random.default_rng
    try:
        ic.np.random.default_rng = lambda: _FixedRNG(np.array([0.0]))
        operator(random_grid, sub_select=mask)
    finally:
        ic.np.random.default_rng = original_default_rng

    assert_allclose(random_grid["density"][:, 0, 0], [0.9, 4.0])

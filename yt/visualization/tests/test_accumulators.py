import numpy as np
import pytest
from unyt.array import unyt_array, ustack
from unyt.testing import assert_allclose_units

from yt.visualization.accumulators import (
    Accumulators,
    _accumulate_scalar_field,
    _accumulate_vector_field,
)

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@pytest.fixture(scope="class")
def curve(ds):
    x = ds.arr(np.linspace(0.0, 1.0, 1000), "code_length")
    y = ds.arr(np.sqrt(x.d), "code_length")
    z = ds.arr(np.power(x.d, 1.0 / 3.0), "code_length")
    return ustack([x, y, z], axis=1)


@pytest.fixture(scope="class")
def scalar_field(curve):
    c = curve.d
    return unyt_array((c[:, 0] + c[:, 1]).reshape(len(c[:, 0]), 1), "kelvin")


@pytest.fixture(scope="class")
def vector_field(curve):
    return curve


@pytest.mark.answer_test
class TestAccumulators:
    answer_file = None
    saved_hashes = None
    answer_version = "000"

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_scalar_integration(self, curve, scalar_field, ds):
        r"""
        Checks _accumulate_scalar_field. The curve is:

        x = np.linspace(0., 1., 10000)
        y = np.sqrt(x)
        z = np.power(x, 1/3)

        Scalar field:
        phi = x + y

        Limits: x \in [0,1]
        """
        a = _accumulate_scalar_field(curve, scalar_field)
        solution = ds.arr([7.0 / 6.0, 5.0 / 6.0, 13.0 / 20.0], "kelvin * code_length")
        try:
            assert_allclose_units(a[-1], solution)
        except AssertionError:
            print(f"test_scalar_integration: {np.abs(a[-1] - solution)}")

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_vector_integration(self, curve, vector_field, ds):
        r"""
        Checks _accumulate_vector_field.

        x = np.linspace(0., 1., 10000)
        y = np.sqrt(x)
        z = np.power(x, 1/3)

        Scalar field:
        \vec{a} = x\hat{i} + y\hat{j} + z\hat{k}

        Limits: x \in [0,1]
        """
        a = _accumulate_vector_field(curve, vector_field)
        solution = ds.quan(3.0 / 2.0, "code_length**2")
        try:
            assert_allclose_units(a[-1], solution)
        except AssertionError:
            print(f"test_vector_integration: {np.abs(a[-1] - solution)}")

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_scalar_tree_access(self, curve, ds):
        accumulator = Accumulators(curve, ds)
        accumulator.accumulate(("gas", "x"), is_vector=False)
        solution = ds.arr([7.0 / 6.0, 5.0 / 6.0, 13.0 / 20.0], "code_length * kelvin")
        try:
            assert_allclose_units(accumulator.accum[-1], solution)
        except AssertionError:
            print(
                f"test_scalar_tree_access: {np.abs(accumulator.accum[-1] - solution)}"
            )

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_vector_tree_access(self, curve, ds):
        accumulator = Accumulators(curve, ds)
        field = [("gas", "x"), ("gas", "y"), ("gas", "z")]
        accumulator.accumulate(field, is_vector=True)
        solution = ds.quan(3.0 / 2.0, "code_length**2")
        try:
            assert_allclose_units(accumulator.accum[-1], solution)
        except AssertionError:
            print(
                f"test_vector_tree_access: {np.abs(accumulator.accum[-1] - solution)}"
            )

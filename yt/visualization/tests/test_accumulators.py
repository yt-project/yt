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


@pytest.fixture(scope="function")
def curve(request, ds):
    """
    The params given to this fixture define the upper limit of
    integration.
    """
    # The upper limit is different for the tree-related and non-tree
    # related tests because, in the tree-related tests, [1., 1., 1.]
    # isn't actually in a node, so there isn't any data
    if "tree" in request.node.name:
        upper_limit = 0.9
    else:
        upper_limit = 1.0
    x = ds.arr(np.linspace(0.0, upper_limit, 10000), "code_length")
    y = ds.arr(np.sqrt(x.d), "code_length")
    z = ds.arr(np.power(x.d, 1.0 / 3.0), "code_length")
    return ustack([x, y, z], axis=1)


@pytest.fixture(scope="function")
def scalar_field(curve):
    c = curve.d
    return unyt_array((c[:, 0] + c[:, 1]).reshape(len(c[:, 0]), 1), "kelvin")


@pytest.fixture(scope="function")
def vector_field(curve):
    return curve


@pytest.fixture(scope="function")
def solution(request, curve, ds):
    # Used in the non-tree related tests
    if curve[-1][0] == 1.0:
        if "scalar" in request.node.name:
            return ds.arr([7.0 / 6.0, 5.0 / 6.0, 13.0 / 20.0], "kelvin * code_length")
        elif "vector" in request.node.name:
            return ds.quan(3.0 / 2.0, "code_length**2")
        else:
            return pytest.fail("Uknown test!")
    # Used in the tree related tests
    elif curve[-1][0] == 0.9:
        if "scalar" in request.node.name:
            return ds.arr(
                [
                    81.0 / 200.0,
                    np.sqrt(729.0 / 1000.0) / 3.0,
                    np.power(6561.0 / 10000.0, 1.0 / 3.0) / 4.0,
                ],
                "code_length**2",
            )
        elif "vector" in request.node.name:
            # Fraction created by the julia function rationalize
            return ds.quan(
                (171.0 / 200.0) + np.power(81.0 / 100.0, 1.0 / 3.0) / 2.0,
                "code_length**2",
            )
        else:
            return pytest.fail("Uknown test!")
    else:
        return pytest.fail("Uknown case!")


@pytest.mark.answer_test
class TestAccumulators:
    answer_file = None
    saved_hashes = None
    answer_version = "000"

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_scalar_integration(self, curve, scalar_field, ds, solution):
        r"""
        Checks _accumulate_scalar_field. The curve is:

        Upper limit is given by the parameters passed to the curve
        fixture
        x = np.linspace(0., upper_limit, 10000)
        y = np.sqrt(x)
        z = np.power(x, 1/3)

        Scalar field:
        phi = x + y

        Limits: x \in [0,1]
        """
        a = _accumulate_scalar_field(curve, scalar_field)
        assert_allclose_units(a[-1], solution, 1e-4)

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_vector_integration(self, curve, vector_field, ds, solution):
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
        assert_allclose_units(a[-1], solution, 1e-4)

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_scalar_tree_access(self, curve, ds, solution):
        accumulator = Accumulators(curve, ds)
        accumulator.accumulate(("gas", "x"), is_vector=False)
        assert_allclose_units(solution, accumulator.accum[0][-1], atol=0.0024)

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_vector_tree_access(self, curve, ds, solution):
        accumulator = Accumulators(curve, ds)
        field = [("gas", "x"), ("gas", "y"), ("gas", "z")]
        accumulator.accumulate(field, is_vector=True)
        assert_allclose_units(solution, accumulator.accum[0][-1], atol=0.0023)

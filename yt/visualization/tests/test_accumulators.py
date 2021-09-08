import numpy as np
import pytest
from unyt.array import unyt_array
from unyt.testing import assert_allclose_units

from yt.visualization.accumulators import (  # Accumulators,
    _accumulate_scalar_field,
    _accumulate_vector_field,
)

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@pytest.fixture(scope="class")
def curve(request):
    x = np.linspace(0.0, 1.0, 10000)
    y = np.sqrt(x)
    z = np.array([0.5] * len(x))
    return unyt_array(np.stack([x, y, z], axis=1), "m")


@pytest.fixture(scope="class")
def scalar_field(curve):
    c = curve.d
    return unyt_array((c[:, 0] + c[:, 1]).reshape(len(c[:, 0]), 1), "kelvin")


@pytest.fixture(scope="class")
def vector_field(curve):
    return unyt_array(curve.d, "m/s")


@pytest.mark.answer_test
class TestAccumulators:
    answer_file = None
    saved_hashes = None
    answer_version = "000"

    def test_scalar(self, curve, scalar_field):
        r"""
        Checks _accumulate_scalar_field. The curve is:

        x = np.linspace(0., 1., 10000)
        y = np.sqrt(x)
        z = 0.5

        Scalar field:
        phi = x + y

        Limits: x \in [0,1]
        """
        a = _accumulate_scalar_field(curve, scalar_field)
        solution = unyt_array([7.0 / 6.0, 5.0 / 6.0, 0.0], "kelvin*m")
        assert_allclose_units(a[-1], solution, 1e-3)

    def test_vector(self, curve, vector_field):
        r"""
        Checks _accumulate_vector_field.

        x = np.linspace(0., 1., 10000)
        y = np.sqrt(x)
        z = 0.5

        Scalar field:
        \vec{a} = x\hat{i} + y\hat{j} + z\hat{k}

        Limits: x \in [0,1]
        """
        a = _accumulate_vector_field(curve, vector_field)
        solution = unyt_array([1.0], "m**2 / s")
        assert_allclose_units(a[-1], solution, 1e-3)

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_same_cell_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in the
        same cell.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_diff_cell_same_node_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in the
        same node but different cells.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_diff_nodes_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in
        different nodes.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_npts_scalar(self, ds):
        r"""
        Calculates the accumulation of the density field along a path with
        more than two points.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_same_cell_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in the same
        cell in the same node.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_diff_cell_same_node_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in the same
        node but in different cells.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_two_pts_diff_nodes_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in different
        nodes.
        """
        raise NotImplementedError

    @pytest.mark.parametrize("ds", [g30], indirect=True)
    def test_npts_vector(self, ds):
        r"""
        Calculates the accumulation of a vector field alog a path with
        more than two points.
        """
        raise NotImplementedError

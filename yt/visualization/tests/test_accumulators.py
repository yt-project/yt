import numpy as np
import pytest

from yt.visualization.accumulators import Accumulators

g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@pytest.fixture(scope="module")
def fake_path(request):
    """
    Samples N (x,y,z) points along y^2=x and z^3=x from (0.1, 0.1, 0.1) ->
    (0.9,0.9^(0.5), 0.9, (1./3.)).
    """
    N = 50
    x = np.linspace(0.1, 0.9, N)
    y = np.linspace(0.1, np.sqrt(0.9), N)
    z = np.linspace(0.1, np.cbrt(0.9), N)
    path = np.stack([x, y, z], axis=1)
    return path


@pytest.mark.answer_test
@pytest.mark.parametrize("ds", [g30], indirect=True)
class TestAccumulators:
    answer_file = None
    saved_hashes = None
    answer_version = "000"

    def test_two_pts_same_cell_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in the same
        cell in the same node.
        """
        path = np.array([[0.16, 0.14, 0.55], [0.155, 0.13, 0.56]])
        accumulator = Accumulators([path], ds)
        field = [("enzo", "Density")]
        accumulator.accumulate(field, is_vector=False)
        answer = np.array([[0.00270012]])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_two_pts_diff_cell_same_node_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in the same
        node but in different cells.
        """
        path = np.array([[0.16, 0.14, 0.44], [0.22, 0.4, 0.8]])
        accumulator = Accumulators([path], ds)
        field = [("enzo", "Density")]
        accumulator.accumulate(field, is_vector=False)
        answer = np.array([[0.080629]])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_two_pts_diff_nodes_scalar(self, ds):
        r"""
        Integrates the density field between two points that are in different
        nodes.
        """
        path = np.array([[0.16, 0.14, 0.55], [0.3, 0.6, 0.75]])
        accumulator = Accumulators([path], ds)
        field = [("enzo", "Density")]
        accumulator.accumulate(field, is_vector=False)
        answer = np.array([[0.09374269]])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_npts_scalar(self, ds, fake_path):
        r"""
        Calculates the accumulation of the density field along the path
        y^2=x and z^3=x from (0.1, 0.1, 0.1) -> (0.9, 0.9^(0.5), 0.9, (1./3.))
        through the IsolatedGalaxy ds.

        NOTE: The answer will depend on the number of points used to define
        the path, so if that changes, update the answer accordingly!
        The current answer is for N = 10
        """
        path = fake_path
        accumulator = Accumulators([path], ds)
        field = [("enzo", "Density")]
        accumulator.accumulate(field, is_vector=False)
        answer = np.array(
            [
                [0.02904837],
                [0.05809688],
                [0.08693819],
                [0.11599762],
                [10.42846079],
                [11.73346686],
                [11.76252416],
                [11.79157306],
                [11.82062148],
            ]
        )
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_two_pts_same_cell_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in the same
        cell in the same node.
        """
        path = np.array([[0.16, 0.14, 0.55], [0.155, 0.13, 0.56]])
        accumulator = Accumulators([path], ds)
        field = [("enzo", "x-velocity"), ("enzo", "y-velocity"), ("enzo", "z-velocity")]
        accumulator.accumulate(field, is_vector=True)
        answer = np.array([-0.00148564])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_two_pts_diff_cell_same_node_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in the same
        node but in different cells.
        """
        path = np.array([[0.16, 0.14, 0.44], [0.22, 0.4, 0.8]])
        field = [("enzo", "x-velocity"), ("enzo", "y-velocity"), ("enzo", "z-velocity")]
        accumulator = Accumulators([path], ds)
        accumulator.accumulate(field, is_vector=True)
        answer = np.array([0.0264723435])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_two_pts_diff_nodes_vector(self, ds):
        r"""
        Integrates the velocity field between two points that are in different
        nodes.
        """
        path = np.array([[0.16, 0.14, 0.55], [0.3, 0.6, 0.75]])
        accumulator = Accumulators([path], ds)
        field = [("enzo", "x-velocity"), ("enzo", "y-velocity"), ("enzo", "z-velocity")]
        accumulator.accumulate(field, is_vector=True)
        answer = np.array([0.05266109])
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

    def test_npts_vector(self, ds, fake_path):
        r"""
        Calculates the accumulation of a test vector field and a test scalar
        field along the path y^2=x and z^3=x from (0.1, 0.1, 0.1) -> (0.9,
        0.9^(0.5), 0.9, (1./3.)) through the IsolatedGalaxy ds.

        NOTE: The answer will depend on the number of points used to define
        the path, so if that changes, update the answer accordingly!
        The current answer is for N = 10
        """
        path = fake_path
        field = [("enzo", "x-velocity"), ("enzo", "y-velocity"), ("enzo", "z-velocity")]
        accumulator = Accumulators([path], ds)
        accumulator.accumulate(field, is_vector=True)
        answer = np.array(
            [
                0.01137986,
                0.03029503,
                0.06936842,
                0.19213309,
                0.96494272,
                0.63585663,
                0.56452914,
                0.5376544,
                0.52296114,
            ]
        )
        np.testing.assert_array_almost_equal(accumulator.accum, answer)

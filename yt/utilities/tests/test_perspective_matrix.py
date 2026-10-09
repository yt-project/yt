import numpy as np
import pytest
from numpy.testing import assert_allclose

from yt.utilities.math_utils import get_perspective_matrix


def _project(matrix, point):
    clip = matrix @ np.append(point, 1.0)
    return clip[:3] / clip[3]


@pytest.mark.parametrize("aspect", [0.5, 1.0, 16 / 9])
def test_perspective_matrix_frustum_corners(aspect):
    fovy, near, far = 60.0, 0.1, 10.0
    matrix = get_perspective_matrix(fovy, aspect, near, far).astype("f8")
    top = near * np.tan(np.radians(fovy) / 2)
    right = top * aspect
    # the corners of the near plane map to the corners of normalized device
    # coordinates, whatever the aspect ratio
    for sx in (-1, 1):
        for sy in (-1, 1):
            ndc = _project(matrix, [sx * right, sy * top, -near])
            assert_allclose(ndc, [sx, sy, -1], atol=1e-6)
    # and the axis stays centered at every depth
    assert_allclose(_project(matrix, [0.0, 0.0, -far])[:2], [0, 0], atol=1e-6)
    assert_allclose(_project(matrix, [0.0, 0.0, -far])[2], 1, atol=1e-5)

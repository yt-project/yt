import numpy as np
from nose.plugins.attrib import attr

from yt.loaders import load_uniform_grid
from yt.testing import ANSWER_TEST_TAG
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import SlicePlot


@attr(ANSWER_TEST_TAG)
def test_inf_and_finite_values_zlim():
    # see https://github.com/yt-project/yt/issues/3901
    shape = (32, 16, 1)
    a = np.ones(16)
    b = np.ones((32, 16))
    c = np.reshape(a * b, shape)

    # injecting an inf
    c[0, 0, 0] = np.inf

    data = {("gas", "density"): c}

    ds = load_uniform_grid(
        data,
        shape,
        bbox=np.array([[0, 1], [0, 1], [0, 1]]),
    )

    def create_image(filename_prefix):
        p = SlicePlot(ds, "z", ("gas", "density"))

        # setting zlim manually
        p.set_zlim(("gas", "density"), -10, 10)
        p.save(filename_prefix)

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_inf_and_finite_values_zlim"
    test.answer_name = "inf_and_finite_values_zlim"
    yield test

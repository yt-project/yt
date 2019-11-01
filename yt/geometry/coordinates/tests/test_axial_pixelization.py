from collections import OrderedDict

import pytest

from yt.testing import \
    fake_amr_ds, _geom_transforms
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAxialPixelization(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    def test_axial_pixelization(self, geom, axis):
        ds = fake_amr_ds(geometry=geom)
        ap = self.axial_pixelization_test(ds)
        if axis == 'x_axis':
            apx_hd = ap[0]
        elif axis == 'y_axis':
            apy_hd = ap[1]
        self.hashes.update({'axial_pixelization' : apx_hd})

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
    def test_axial_pixelization(self):
        self.hashes['axial_pixelization'] = OrderedDict()
        self.hashes['axial_pixelization']['x_axis'] = OrderedDict()
        self.hashes['axial_pixelization']['y_axis'] = OrderedDict()
        for geom in sorted(_geom_transforms):
            ds = fake_amr_ds(geometry=geom)
            ap = self.axial_pixelization_test(ds)
            apx_hd = ap[0]
            apy_hd = ap[1]
            self.hashes['axial_pixelization']['x_axis'][geom] = apx_hd
            self.hashes['axial_pixelization']['y_axis'][geom] = apy_hd

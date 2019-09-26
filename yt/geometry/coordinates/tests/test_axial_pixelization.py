from collections import OrderedDict

import pytest

from yt.testing import \
    fake_amr_ds, _geom_transforms
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestAxialPixelization(fw.AnswerTest):
    def test_axial_pixelization(self):
        hd = OrderedDict()
        hd['axial_pixelization'] = OrderedDict()
        hd['axial_pixelization']['x_axis'] = OrderedDict()
        hd['axial_pixelization']['y_axis'] = OrderedDict()
        for geom in sorted(_geom_transforms):
            ds = fake_amr_ds(geometry=geom)
            ap = self.axial_pixelization_test(ds)
            apx_hd = utils.generate_hash(ap[0])
            apy_hd = utils.generate_hash(ap[1])
            hd['axial_pixelization']['x_axis'][geom] = apx_hd
            hd['axial_pixelization']['y_axis'][geom] = apy_hd
        hashes = {'axial_pixelization' : hd}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

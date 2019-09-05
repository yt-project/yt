import pytest

from yt.testing import \
    fake_amr_ds, _geom_transforms
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestAxialPixelization(fw.AnswerTest):
    def test_axial_pixelization(self):
        apx_hd = b''
        apy_hd = b''
        for geom in sorted(_geom_transforms):
            ds = fake_amr_ds(geometry=geom)
            ap = self.axial_pixelization_test(ds)
            apx_hd += ap[0]
            apy_hd += ap[1]
        hashes = {'axial_pixelization_x' : utils.generate_hash(apx_hd),
            'axial_pixelization_y' : utils.generate_hash(apy_hd)}
        utils.handle_hashes(self.save_dir, 'axial_pixelization', hashes, self.answer_store)

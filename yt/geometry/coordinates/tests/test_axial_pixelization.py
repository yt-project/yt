import pytest

from yt.testing import fake_amr_ds 
import yt.utilities.answer_testing.framework as fw


@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestAxialPixelization(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    def test_axial_pixelization(self, geom, axis):
        ds = fake_amr_ds(geometry=geom)
        ap = self.axial_pixelization_test(ds)
        if axis == 'x_axis':
            ap_hd = ap[0]
        elif axis == 'y_axis':
            ap_hd = ap[1]
        self.hashes.update({'axial_pixelization' : ap_hd})

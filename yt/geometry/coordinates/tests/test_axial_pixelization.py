import pytest

from yt.testing import fake_amr_ds
from yt.utilities.answer_testing.answer_tests import axial_pixelization


@pytest.mark.answer_test
class TestAxialPixelization:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    def test_axial_pixelization(self, geom):
        ds = fake_amr_ds(geometry=geom)
        ap = axial_pixelization(ds)
        self.hashes.update({"axial_pixelization": ap})

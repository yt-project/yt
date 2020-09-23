import pytest

import yt
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import generic_image


def compare(ds, field):
    def slice_image(im_name):
        sl = yt.SlicePlot(ds, "z", field)
        sl.set_log("all", False)
        image_file = sl.save(im_name)
        return image_file

    gi = generic_image(slice_image)
    # generic_image returns a list. In this case, there's only one entry,
    # which is a np array with the data we want
    return gi[0]


raw_fields = "Laser/plt00015"


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir")
class TestRawFieldSlices:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    @utils.requires_ds(raw_fields)
    def test_raw_field_slices(self, field):
        ds = utils.data_dir_load(raw_fields)
        gi = compare(ds, field)
        self.hashes.update({"generic_image": gi})

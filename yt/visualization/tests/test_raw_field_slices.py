"""
Tests for making slices through raw fields

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
import pytest

import yt
from yt.utilities.answer_testing import utils
from yt.utilities.answer_testing.answer_tests import generic_image

raw_fields = "Laser/plt00015"


@pytest.mark.answer_test
@pytest.mark.usefixtures("temp_dir")
class TestRawFieldSlices:
    @pytest.mark.usefixtures("hashing")
    @utils.requires_ds(raw_fields)
    def test_raw_field_slices(self, field):
        ds = utils.data_dir_load(raw_fields)
        sl = yt.SlicePlot(ds, "z", field)
        sl.set_log("all", False)
        image_file = sl.save("slice_answers_raw_{}".format(field[1]))
        gi = generic_image(image_file[0])
        self.hashes.update({"generic_image": gi})

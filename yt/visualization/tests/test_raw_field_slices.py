"""
Tests for making slices through raw fields

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from collections import OrderedDict

import pytest

import yt
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

raw_fields = "Laser/plt00015"
_raw_field_names =  [('raw', 'Bx'),
                     ('raw', 'By'),
                     ('raw', 'Bz'),
                     ('raw', 'Ex'),
                     ('raw', 'Ey'),
                     ('raw', 'Ez'),
                     ('raw', 'jx'),
                     ('raw', 'jy'),
                     ('raw', 'jz')]

@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestRawFieldSlices(fw.AnswerTest):
    @utils.requires_ds(raw_fields)
    def test_raw_field_slices(self):
        ds = utils.data_dir_load(raw_fields)
        hd = OrderedDict()
        hd['generic_image'] = OrderedDict()
        for field in _raw_field_names:
            sl = yt.SlicePlot(ds, 'z', field)
            sl.set_log('all', False)
            image_file = sl.save("slice_answers_raw_{}".format(field[1]))
            gi_hd = utils.generate_hash(self.generic_image_test(image_file))
            hd['generic_image'][field] = gi_hd 
        hashes = {'raw-field-slices' : hd} 
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store) 

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

import yt
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
@pytest.mark.usefixtures('temp_dir')
@requires_ds(raw_fields)
def test_raw_field_slices():
    ds = utils.data_dir_load(raw_fields)
    for field in _raw_field_names:
        sl = yt.SlicePlot(ds, 'z', field)
        sl.set_log('all', False)
        image_file = sl.save("slice_answers_raw_{}".format(field[1]))
        gi_hd += self.generic_image_test(image_file)
    hashes = {'raw-field-slices' : utils.generate_hash(gi_hd)}
    utils.handle_hashes(self.save_dir, 'raw-field_slices', hashes, self.answer_store) 

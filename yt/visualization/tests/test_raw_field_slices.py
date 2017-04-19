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
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    data_dir_load, \
    GenericImageTest


def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def compare(ds, field, test_prefix, decimals=12):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, 'z', field)
        sl.set_log('all', False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = "slice_{}".format(test_prefix)
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    return test


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
@requires_ds(raw_fields)
def test_raw_field_slices():
    ds = data_dir_load(raw_fields)
    for field in _raw_field_names:
        yield compare(ds, field, "answers_raw_%s" % field[1])


"""
Arepo frontend tests




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import OrderedDict
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load, \
    requires_ds, \
    sph_answer
from yt.frontends.arepo.api import ArepoHDF5Dataset

bullet_h5 = "ArepoBullet/snapshot_150.hdf5"

# py2/py3 compat
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

@requires_file(bullet_h5)
def test_arepo_hdf5():
    assert isinstance(data_dir_load(bullet_h5),
                      ArepoHDF5Dataset)


bullet_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None)
    ]
)


@requires_ds(bullet_h5)
def test_arepo_bullet():
    ds = data_dir_load(bullet_h5)
    for test in sph_answer(ds, 'snapshot_150', 40193369, 
                           bullet_fields):
        test_arepo_bullet.__name__ = test.description
        yield test

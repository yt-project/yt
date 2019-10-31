"""
HaloCatalog answer tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os

import numpy as np
import pytest

from yt.analysis_modules.halo_analysis.api import \
    HaloCatalog, \
    add_quantity
from yt.convenience import \
    load
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

rh0 = "rockstar_halos/halos_0.0.bin"
e64 = "Enzo_64/DD0043/data0043"


def _nstars(halo):
    sp = halo.data_object
    return (sp["all", "creation_time"] > 0).sum()
add_quantity("nstars", _nstars)

@pytest.mark.answer_test
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestHaloQuantity(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(rh0)
    @utils.requires_ds(e64)
    def test_halo_quantity(self):
        data_ds_fn = e64
        halos_ds_fn = rh0
        ds = utils.data_dir_load(data_ds_fn)

        dds = utils.data_dir_load(data_ds_fn)
        hds = utils.data_dir_load(halos_ds_fn)
        hc = HaloCatalog(
            data_ds=dds, halos_ds=hds,
            output_dir=os.path.join(os.getcwd(), str(dds)))
        hc.add_callback("sphere")
        hc.add_quantity("nstars")
        hc.create()

        fn = os.path.join(os.getcwd(), str(dds),
                          "%s.0.h5" % str(dds))
        ds = load(fn)
        ad = ds.all_data()
        mi, ma = ad.quantities.extrema("nstars")
        mean = ad.quantities.weighted_average_quantity(
            "nstars", "particle_ones")
        self.hashes['halo_quantity'] = np.array([mean, mi, ma])

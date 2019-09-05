"""
light cone generator test



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import pytest

from yt.testing import \
    requires_module
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


ETC = "enzo_tiny_cosmology/32Mpc_32.enzo"


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir')
@requires_module("h5py")
@utils.requires_sim(ETC, "Enzo")
def test_light_cone_projection():
    lcp_hd = self.light_cone_projection_test(ETC, "Enzo")
    hashes = {'light_cone_projection' : utils.generate_hash(lcp_hd)}
    utils.handle_hashes(self.save_dir, 'light_cone', hashes, self.answer_store)

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


@pytest.mark.answer_test
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestLightCone(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @requires_module("h5py")
    @utils.requires_sim(ETC, "Enzo")
    def test_light_cone_projection(self):
        lcp_hd = self.light_cone_projection_test(ETC, "Enzo")
        self.hashes.update({'light_cone_projection' : lcp_hd})

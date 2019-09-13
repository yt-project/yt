"""
Unit test for the RADMC3D Exporter analysis module
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import shutil
import tempfile

import numpy as np
import pytest

import yt
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter
from yt.testing import assert_allclose
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


ISO_GAL = "IsolatedGalaxy/galaxy0030/galaxy0030"


# Answer file
answer_file = 'radmc3d_answers.yaml'


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
reason="--with-answer-testing not set.")
class TestRadmc3dExporter(fw.AnswerTest):
    @pytest.mark.usefixtures('temp_dir')
    @utils.requires_ds(ISO_GAL)
    def test_radmc3d_exporter_continuum(self):
        """
        This test is simply following the description in the docs for how to
        generate the necessary output files to run a continuum emission map from
        dust for one of our sample datasets.
        """
        ds = yt.load(ISO_GAL)
        # Make up a dust density field where dust density is 1% of gas density
        dust_to_gas = 0.01
        def _DustDensity(field, data):
            return dust_to_gas * data["density"]
        ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")
        # try to write the output files
        writer = RadMC3DWriter(self.ds)
        writer.write_amr_grid()
        writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
        # compute the sum of the values in the resulting file
        total = 0.0
        with open('dust_density.inp', 'r') as f:
            for i, line in enumerate(f):
                # skip header
                if i < 3:
                    continue
                line = line.rstrip()
                total += np.float64(line)
        total = total.tostring()
        hashes = {'radmc3d_exporter_continuum' : utils.generate_hash(total)}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

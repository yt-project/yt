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

import yt
from yt.testing import assert_allclose
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter
from yt.utilities.answer_testing.framework import \
    AnswerTestingTest, \
    requires_ds
import tempfile
import numpy as np
import os
import shutil


class RadMC3DValuesTest(AnswerTestingTest):
    '''

    This test writes out a "dust_density.inp" file, 
    reads it back in, and checks the sum of the 
    values for degradation.

    '''
    _type_name = "RadMC3DValuesTest"
    _attrs = ("field", )

    def __init__(self, ds_fn, field, decimals=10):
        super(RadMC3DValuesTest, self).__init__(ds_fn)
        self.field = field
        self.decimals = decimals

    def run(self):

        # Set up in a temp dir
        tmpdir = tempfile.mkdtemp()
        curdir = os.getcwd()
        os.chdir(tmpdir)

        # try to write the output files
        writer = RadMC3DWriter(self.ds)
        writer.write_amr_grid()
        writer.write_dust_file(self.field, "dust_density.inp")

        # compute the sum of the values in the resulting file
        total = 0.0
        with open('dust_density.inp', 'r') as f:
            for i, line in enumerate(f):

                # skip header
                if i < 3:
                    continue

                line = line.rstrip()
                total += np.float64(line)

        # clean up
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

        return total

    def compare(self, new_result, old_result):
        err_msg = "Total value for %s not equal." % (self.field,)
        assert_allclose(new_result, old_result, 10.**(-self.decimals),
                        err_msg=err_msg, verbose=True)


ISO_GAL = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_ds(ISO_GAL)
def test_radmc3d_exporter_continuum():
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

    yield RadMC3DValuesTest(ds, ("gas", "dust_density"))

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
from yt.testing import *
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter
from yt.config import ytcfg
import tempfile
import os
import shutil

ISO_GAL = "IsolatedGalaxy/galaxy0030/galaxy0030"

@requires_file(ISO_GAL)
def test_radmc3d_exporter_continuum():
    """
    This test is simply following the description in the docs for how to
    generate the necessary output files to run a continuum emission map from
    dust for one of our sample datasets.
    """

    # Set up in a temp dir
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # Make up a dust density field where dust density is 1% of gas density
    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data["density"]
    
    # Make up a dust temperature field where temp = 10K everywhere
    def _DustTemperature(field, data):
        return 0.*data["temperature"] + data.ds.quan(10,'K')
    
    ds = yt.load(ISO_GAL)
    ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")
    ds.add_field(("gas", "dust_temperature"), function=_DustTemperature, units="K")
    writer = RadMC3DWriter(ds)
    
    writer.write_amr_grid()
    writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
    writer.write_dust_file(("gas", "dust_temperature"), "dust_temperature.inp")

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

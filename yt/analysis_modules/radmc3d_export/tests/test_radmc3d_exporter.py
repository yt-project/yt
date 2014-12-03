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
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter
from yt.units.yt_array import YTQuantity
from yt.testing import *

def test_radmc3d_exporter():
    """
    This test is simply following the description in the docs for how to
    generate the necessary output files to run a continuum emission map from
    dust for one of our sample datasets.
    """

    # Make up a dust density field where dust density is 1% of gas density
    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data["density"]
    yt.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")
    
    # Make up a dust temperature field where temp = 10K everywhere
    def _DustTemperature(field, data):
        return 0.*data["temperature"] + data.ds.quan(10,'K')
    yt.add_field(("gas", "dust_temperature"), function=_DustTemperature, units="K")
    
    ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    writer = RadMC3DWriter(ds)
    
    writer.write_amr_grid()
    writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
    writer.write_dust_file(("gas", "dust_temperature"), "dust_temperature.inp")

"""
Answer test the photon_simulator analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import *
from yt.analysis_modules.photon_simulator.api import *
from yt.utilities.answer_testing.framework import requires_pf, \
     GenericArrayTest, data_dir_load
import numpy as np

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

MHD = "MHDSloshing/virgo_low_res.0054.vtk"
@requires_module("xspec")
@requires_pf(MHD)
def test_cluster():
    np.random.seed(seed=0x4d3d3d3)
    pf = data_dir_load(MHD, parameters={"TimeUnits":3.1557e13,
                                        "LengthUnits":3.0856e24,
                                        "DensityUnits":6.770424595218825e-27})

    A = 3000.
    exp_time = 3.0e5
    redshift = 0.02
    
    apec_model = XSpecThermalModel("apec", 0.1, 20.0, 2000)
    tbabs_model = XSpecAbsorbModel("TBabs", 0.1)

    ARF = os.environ["YT_DATA_DIR"]+"/xray_data/chandra_ACIS-S3_onaxis_arf.fits"
    RMF = os.environ["YT_DATA_DIR"]+"/xray_data/chandra_ACIS-S3_onaxis_rmf.fits"
            
    sphere = pf.h.sphere("c", (0.25, "mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model, cosmology=cosmo)
    
    events = photons.project_photons([0.0,0.0,1.0],
                                     responses=[ARF,RMF],
                                     absorb_model=abs_model)
    
    for k,v in photons.items():
        def photons_test(v): return v
        test = GenericArrayTest(pf, photons_test)
        test_cluster.__name__ = test.description
        yield test

    for k,v in events.items():
        def events_test(v): return v
        test = GenericArrayTest(pf, events_test)
        test_cluster.__name__ = test.description
        yield test
            
            

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
from yt.config import ytcfg
from yt.analysis_modules.photon_simulator.api import *
from yt.utilities.answer_testing.framework import requires_ds, \
     GenericArrayTest, data_dir_load
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

test_dir = ytcfg.get("yt", "test_data_dir")

ETC = test_dir+"/enzo_tiny_cosmology/DD0046/DD0046"
APEC = test_dir+"/xray_data/atomdb_v2.0.2"
TBABS = test_dir+"/xray_data/tbabs_table.h5"
ARF = test_dir+"/xray_data/chandra_ACIS-S3_onaxis_arf.fits"
RMF = test_dir+"/xray_data/chandra_ACIS-S3_onaxis_rmf.fits"

@requires_ds(ETC)
@requires_file(APEC)
@requires_file(TBABS)
@requires_file(ARF)
@requires_file(RMF)
def test_etc():

    np.random.seed(seed=0x4d3d3d3)

    ds = data_dir_load(ETC)
    A = 3000.
    exp_time = 1.0e5
    redshift = 0.1

    apec_model = TableApecModel(APEC, 0.1, 20.0, 2000)
    tbabs_model = TableAbsorbModel(TBABS, 0.1)

    sphere = ds.sphere("max", (0.5, "mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    events = photons.project_photons([0.0,0.0,1.0],
                                     responses=[ARF,RMF],
                                     absorb_model=tbabs_model)

    def photons_test(): return photons.photons
    def events_test(): return events.events

    for test in [GenericArrayTest(ds, photons_test),
                 GenericArrayTest(ds, events_test)]:
        test_etc.__name__ = test.description
        yield test

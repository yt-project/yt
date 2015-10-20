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

from yt.analysis_modules.photon_simulator.api import \
    TableApecModel, TableAbsorbModel, \
    ThermalPhotonModel, PhotonList
from yt.config import ytcfg
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import requires_ds, \
    GenericArrayTest, data_dir_load
import numpy as np

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

test_data_dir = ytcfg.get("yt", "test_data_dir")
xray_data_dir = ytcfg.get("yt", "xray_data_dir")

ETC = test_data_dir+"/enzo_tiny_cosmology/DD0046/DD0046"
APEC = xray_data_dir+"/atomdb_v2.0.2"
TBABS = xray_data_dir+"/tbabs_table.h5"
ARF = xray_data_dir+"/aciss_aimpt_cy17.arf"
RMF = xray_data_dir+"/aciss_aimpt_cy17.rmf"

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

    sphere = ds.sphere("max", (0.5, "Mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    events = photons.project_photons("z", responses=[ARF,RMF],
                                     absorb_model=tbabs_model)

    def photons_test():
        return photons.photons
    def events_test():
        return events.events

    for test in [GenericArrayTest(ds, photons_test),
                 GenericArrayTest(ds, events_test)]:
        test_etc.__name__ = test.description
        yield test

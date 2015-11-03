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
from numpy.random import RandomState
import os

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

test_data_dir = ytcfg.get("yt", "test_data_dir")
xray_data_dir = ytcfg.get("yt", "xray_data_dir")

rmfs = ["pn-med.rmf", "acisi_aimpt_cy17.rmf",
        "aciss_aimpt_cy17.rmf", "nustar.rmf",
        "ah_sxs_5ev_basefilt_20100712.rmf"]
arfs = ["pn-med.arf", "acisi_aimpt_cy17.arf",
        "aciss_aimpt_cy17.arf", "nustar_3arcminA.arf",
        "sxt-s_120210_ts02um_intallpxl.arf"]

gslr = test_data_dir+"/GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
APEC = xray_data_dir
TBABS = xray_data_dir+"/tbabs_table.h5"

def return_data(data):
    def _return_data():
        return data
    return _return_data

@requires_ds(gslr)
@requires_file(APEC)
@requires_file(TBABS)
def test_sloshing():

    prng = RandomState(0x4d3d3d3)

    ds = data_dir_load(gslr)
    A = 2000.
    exp_time = 1.0e4
    redshift = 0.1

    apec_model = TableApecModel(APEC, 0.1, 11.0, 10000)
    tbabs_model = TableAbsorbModel(TBABS, 0.1)

    sphere = ds.sphere("c", (0.1, "Mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3, prng=prng)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    return_photons = return_data(photons.photons)

    tests = []
    tests.append(GenericArrayTest(ds, return_photons))

    for a, r in zip(arfs, rmfs):
        arf = os.path.join(xray_data_dir,a)
        rmf = os.path.join(xray_data_dir,r)
        events = photons.project_photons([1.0,-0.5,0.2], responses=[arf,rmf],
                                         absorb_model=tbabs_model, 
                                         convolve_energies=True, prng=prng)

        return_events = return_data(events.events)

        tests.append(GenericArrayTest(ds, return_events))

    for test in tests:
        test_sloshing.__name__ = test.description
        yield test

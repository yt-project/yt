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
    ThermalPhotonModel, PhotonList, EventList, \
    convert_old_file, merge_files
from yt.config import ytcfg
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import requires_ds, \
    GenericArrayTest, data_dir_load
import numpy as np
from numpy.testing import assert_array_equal
from numpy.random import RandomState
import os
import tempfile
import shutil

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

gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
APEC = xray_data_dir
TBABS = os.path.join(xray_data_dir, "tbabs_table.h5")
old_photons = os.path.join(xray_data_dir, "old_photons.h5")
old_events = os.path.join(xray_data_dir, "old_events.h5")

def return_data(data):
    def _return_data(name):
        return data
    return _return_data

@requires_ds(gslr)
@requires_file(APEC)
@requires_file(TBABS)
@requires_file(old_photons)
@requires_file(old_events)
def test_sloshing():

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

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

    tests = [GenericArrayTest(ds, return_photons, args=["photons"])]

    for a, r in zip(arfs, rmfs):
        arf = os.path.join(xray_data_dir, a)
        rmf = os.path.join(xray_data_dir, r)
        events = photons.project_photons([1.0,-0.5,0.2], responses=[arf,rmf],
                                         absorb_model=tbabs_model, 
                                         convolve_energies=True, prng=prng)

        return_events = return_data(events.events)

        tests.append(GenericArrayTest(ds, return_events, args=[a]))

    for test in tests:
        test_sloshing.__name__ = test.description
        yield test

    photons.write_h5_file("test_photons.h5")
    events.write_h5_file("test_events.h5")

    photons2 = PhotonList.from_file("test_photons.h5")
    events2 = EventList.from_h5_file("test_events.h5")

    convert_old_file("old_photons.h5", "converted_photons.h5")
    convert_old_file("old_events.h5", "converted_events.h5")

    photons3 = PhotonList.from_file("converted_photons.h5")
    events3 = EventList.from_h5_file("converted_events.h5")

    for k in photons:
        yield assert_array_equal, photons[k].d, photons2[k].d
        yield assert_array_equal, photons[k].d, photons3[k].d
    for k in events:
        yield assert_array_equal, events[k].d, events2[k].d
        yield assert_array_equal, events[k].d, events3[k].d

    nevents = 0

    for i in range(4):
        events = photons.project_photons([1.0,-0.5,0.2],
                                         exp_time=0.25*exp_time,
                                         absorb_model=tbabs_model,
                                         prng=prng)
        events.write_h5_file("split_events_%d.h5" % i)
        nevents += len(events["xsky"])

    merge_files(["split_events_%d.h5" % i for i in range(4)],
                "merged_events.h5", add_exposure_times=True)

    merged_events = EventList.from_h5_file("merged_events.h5")
    assert len(merged_events["xsky"]) == nevents
    assert merged_events.parameters["ExposureTime"] == exp_time

    os.chdir(curdir)
    shutil.rmtree(tmpdir)

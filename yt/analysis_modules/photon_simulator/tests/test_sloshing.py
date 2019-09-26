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
from collections import OrderedDict
import os

from numpy.random import RandomState
import pytest

from yt.analysis_modules.photon_simulator.api import \
    TableApecModel, TableAbsorbModel, \
    ThermalPhotonModel, PhotonList, EventList, \
    convert_old_file, merge_files
from yt.config import ytcfg
from yt.testing import \
    requires_file, \
    assert_almost_equal
from yt.units.yt_array import uconcatenate
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

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
old_photon_file = os.path.join(xray_data_dir, "old_photons.h5")
old_event_file = os.path.join(xray_data_dir, "old_events.h5")


def return_data(data):
    def _return_data(name):
        return data
    return _return_data


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir')
class TestSloshingPhoton(fw.AnswerTest):
    @utils.requires_ds(gslr)
    @requires_file(APEC)
    @requires_file(TBABS)
    @requires_file(old_photon_file)
    @requires_file(old_event_file)
    def test_sloshing(self):
        prng = RandomState(0x4d3d3d3)
        ds = utils.data_dir_load(gslr)
        A = 2000.
        exp_time = 1.0e4
        redshift = 0.1
        apec_model = TableApecModel(APEC, 0.1, 11.0, 10000)
        tbabs_model = TableAbsorbModel(TBABS, 0.1)
        sphere = ds.sphere("c", (0.1, "Mpc"))
        thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3, prng=prng)
        photons1 = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                           thermal_model)
        return_photons = return_data(photons1.photons)
        hd = OrderedDict()
        ga_hd = utils.generate_hash(
            self.generic_array_test(ds, return_photons, args=['photons'])
        )
        hd['generic_array1'] = ga_hd
        hd['generic_array2'] = OrderedDict()
        for a, r in zip(arfs, rmfs):
            arf = os.path.join(xray_data_dir, a)
            rmf = os.path.join(xray_data_dir, r)
            events1 = photons1.project_photons([1.0,-0.5,0.2], responses=[arf,rmf],
                                              absorb_model=tbabs_model, 
                                              convolve_energies=True, prng=prng)
            events1['xsky']
            return_events = return_data(events1.events)
            ga_hd = utils.generate_hash(
                self.generic_array_test(ds, return_events, args=[a])
            )
            hd['generic_array2'][a] = ga_hd
        hashes = {'sloshing' : hd}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)
        photons1.write_h5_file("test_photons.h5")
        events1.write_h5_file("test_events.h5")
        photons2 = PhotonList.from_file("test_photons.h5")
        events2 = EventList.from_h5_file("test_events.h5")
        convert_old_file(old_photon_file, "converted_photons.h5")
        convert_old_file(old_event_file, "converted_events.h5")
        PhotonList.from_file("converted_photons.h5")
        EventList.from_h5_file("converted_events.h5")
        for k in photons1.keys():
            if k == "Energy":
                arr1 = uconcatenate(photons1[k])
                arr2 = uconcatenate(photons2[k])
            else:
                arr1 = photons1[k]
                arr2 = photons2[k]
            assert_almost_equal(arr1, arr2)
        for k in events1.keys():
            assert_almost_equal(events1[k], events2[k])
        nevents = 0
        for i in range(4):
            events = photons1.project_photons([1.0,-0.5,0.2],
                                             exp_time_new=0.25*exp_time,
                                             absorb_model=tbabs_model,
                                             prng=prng)
            events.write_h5_file("split_events_%d.h5" % i)
            nevents += len(events["xsky"])
        merge_files(["split_events_%d.h5" % i for i in range(4)],
                    "merged_events.h5", add_exposure_times=True,
                    clobber=True)
        merged_events = EventList.from_h5_file("merged_events.h5")
        assert len(merged_events["xsky"]) == nevents
        assert merged_events.parameters["ExposureTime"] == exp_time

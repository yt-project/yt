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
import os

from numpy.random import RandomState
import pytest

from yt.analysis_modules.photon_simulator.api import \
    PhotonList, EventList, TableApecModel, TableAbsorbModel, \
    ThermalPhotonModel, PhotonList, convert_old_file, merge_files
from yt.config import ytcfg
from yt.testing import assert_almost_equal
from yt.units.yt_array import uconcatenate
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


xray_data_dir = ytcfg.get("yt", "xray_data_dir")
gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
old_photon_file = os.path.join(xray_data_dir, "old_photons.h5")
old_event_file = os.path.join(xray_data_dir, "old_events.h5")
TBABS = os.path.join(xray_data_dir, "tbabs_table.h5")


def return_data(data):
    def _return_data(name):
        return data
    return _return_data



@pytest.mark.answer_test
@pytest.mark.usefixtures('temp_dir', 'answer_file', 'hashing')
class TestSloshingPhoton(fw.AnswerTest):
    def test_sloshing_return_photons(self):
        ds = utils.data_dir_load(gslr)
        prng = RandomState(0x4d3d3d3)
        A = 2000.
        exp_time = 1.0e4
        redshift = 0.1
        apec_model = TableApecModel(xray_data_dir, 0.1, 11.0, 10000)
        tbabs_model = TableAbsorbModel(TBABS, 0.1)
        sphere = ds.sphere("c", (0.1, "Mpc"))
        thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3, prng=prng)
        photons1 = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                           thermal_model)
        rdata = return_data(photons1.photons)
        ga_hd = self.generic_array_test(rdata, args=['photons'])
        self.hashes.update({'generic_array' : ga_hd})

    def test_sloshing_return_events(self, arf, rmf):
        arf = os.path.join(xray_data_dir, arf)
        rmf = os.path.join(xray_data_dir, rmf)
        ds = utils.data_dir_load(gslr)
        prng = RandomState(0x4d3d3d3)
        A = 2000.
        exp_time = 1.0e4
        redshift = 0.1
        apec_model = TableApecModel(xray_data_dir, 0.1, 11.0, 10000)
        tbabs_model = TableAbsorbModel(TBABS, 0.1)
        sphere = ds.sphere("c", (0.1, "Mpc"))
        thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3, prng=prng)
        photons1 = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                           thermal_model)
        rdata = return_data(photons1.photons)
        events1 = photons1.project_photons([1.0,-0.5,0.2], responses=[arf,rmf],
          absorb_model=tbabs_model, 
          convolve_energies=True, prng=prng
        )
        events1['xsky']
        return_events = return_data(events1.events)
        ga_hd = self.generic_array_test(return_events, args=[os.path.basename(arf)])
        self.hashes.update({'generic_array_return_events' : ga_hd})
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

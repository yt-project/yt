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

import pytest

from yt.analysis_modules.photon_simulator.api import \
    ThermalPhotonModel, PhotonList, EventList, \
    convert_old_file, merge_files
from yt.testing import \
    requires_file, \
    assert_almost_equal
from yt.units.yt_array import uconcatenate
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
old_photon_file = os.path.join(xray_data_dir, "old_photons.h5")
old_event_file = os.path.join(xray_data_dir, "old_events.h5")


def return_data(data):
    def _return_data(name):
        return data
    return _return_data


@pytest.mark.answer_test
@pytest.mark.usefixtures('temp_dir', 'answer_file')
class TestSloshingPhoton(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    def test_sloshing_return_photons(self, photon_data, ds_gslr):
        ga_hd = self.generic_array_test(ds_gslr, photon_data.return_photons, args=['photons'])
        self.hashes.update({'generic_array' : ga_hd})

    @pytest.mark.usefixtures('hashing')
    def test_sloshing_return_events(self, arf, rmf, photon_data, ds_gslr)
        events1 = photon_data.photons1.project_photons([1.0,-0.5,0.2], responses=[arf,rmf],
          absorb_model=photon_data.tbabs_model, 
          convolve_energies=True, prng=photon_data.prng
        )
        events1['xsky']
        return_events = photon_data.return_data(events1.events)
        ga_hd = self.generic_array_test(ds_gslr, return_events, args=[os.path.basename(arf)])
        self.hashes.upate({'generic_array_return_events' : ga_hd})
        photon_data.photons1.write_h5_file("test_photons.h5")
        events1.write_h5_file("test_events.h5")
        photons2 = PhotonList.from_file("test_photons.h5")
        events2 = EventList.from_h5_file("test_events.h5")
        convert_old_file(old_photon_file, "converted_photons.h5")
        convert_old_file(old_event_file, "converted_events.h5")
        PhotonList.from_file("converted_photons.h5")
        EventList.from_h5_file("converted_events.h5")
        for k in photons1.keys():
            if k == "Energy":
                arr1 = uconcatenate(photon_data.photons1[k])
                arr2 = uconcatenate(photons2[k])
            else:
                arr1 = photon_data.photons1[k]
                arr2 = photons2[k]
            assert_almost_equal(arr1, arr2)
        for k in events1.keys():
            assert_almost_equal(events1[k], events2[k])
        nevents = 0
        for i in range(4):
            events = photon_data.photons1.project_photons([1.0,-0.5,0.2],
                                             exp_time_new=0.25*photon_data.exp_time,
                                             absorb_model=photon_data.tbabs_model,
                                             prng=photon_data.prng)
            events.write_h5_file("split_events_%d.h5" % i)
            nevents += len(events["xsky"])
        merge_files(["split_events_%d.h5" % i for i in range(4)],
                    "merged_events.h5", add_exposure_times=True,
                    clobber=True)
        merged_events = EventList.from_h5_file("merged_events.h5")
        assert len(merged_events["xsky"]) == nevents
        assert merged_events.parameters["ExposureTime"] == exp_time

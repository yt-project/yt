import os

from numpy.random import RandomState
import pytest

from yt.analysis_modules.photon_simulator.api import \
    TableApecModel, TableAbsorbModel, \
    ThermalPhotonModel, PhotonList, EventList, \
    convert_old_file, merge_files
from yt.config import ytcfg
from yt.utilities.answer_testing import utils


gslr = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
test_data_dir = ytcfg.get("yt", "test_data_dir")
xray_data_dir = ytcfg.get("yt", "xray_data_dir")
rmfs = ["pn-med.rmf", "acisi_aimpt_cy17.rmf",
        "aciss_aimpt_cy17.rmf", "nustar.rmf",
        "ah_sxs_5ev_basefilt_20100712.rmf"]
arfs = ["pn-med.arf", "acisi_aimpt_cy17.arf",
        "aciss_aimpt_cy17.arf", "nustar_3arcminA.arf",
        "sxt-s_120210_ts02um_intallpxl.arf"]
fs_ids = ['pn-med', 'acisi_aimpt_cy17', 'aciss_aimpty_cy17', 'nustar_3arcminA',
    'sxt-s_120210_ts02um_intallpxl']
TBABS = os.path.join(xray_data_dir, "tbabs_table.h5")


def return_data(data):
    def _return_data(name):
        return data
    return _return_data


class PhotonData:
    def __init__(self, ds):
        self.prng = RandomState(0x4d3d3d3)
        A = 2000.
        self.exp_time = 1.0e4
        redshift = 0.1
        self.apec_model = TableApecModel(xray_data_dir, 0.1, 11.0, 10000)
        self.tbabs_model = TableAbsorbModel(TBABS, 0.1)
        sphere = ds.sphere("c", (0.1, "Mpc"))
        self.thermal_model = ThermalPhotonModel(self.apec_model, Zmet=0.3, prng=self.prng)
        self.photons1 = PhotonList.from_scratch(sphere, redshift, A, self.exp_time,
                                           self.thermal_model)
        self.return_data = return_data(self.photons1.photons)

@pytest.fixture(scope='function')
def photon_data(request):
    if request.function.__name__ in ['test_sloshing_return_photons', 'test_sloshing_return_events']:
        ds = utils.data_dir_load(gslr)
        return PhotonData(ds)

def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_sloshing_return_events':
        params = [[a, r] for a, r in zip(arfs, rmfs)]
        for i in range(len(params)):
            params[i][0] = os.path.join(xray_data_dir, params[i][0])
            params[i][1] = os.path.join(xray_data_dir, params[i][1])
        metafunc.parametrize('arf, rmf', params, ids=fs_ids)

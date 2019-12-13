import os

import pytest

from yt.config import ytcfg


xray_data_dir = ytcfg.get("yt", "xray_data_dir")


rmfs = ["pn-med.rmf", "acisi_aimpt_cy17.rmf",
        "aciss_aimpt_cy17.rmf", "nustar.rmf",
        "ah_sxs_5ev_basefilt_20100712.rmf"]
arfs = ["pn-med.arf", "acisi_aimpt_cy17.arf",
        "aciss_aimpt_cy17.arf", "nustar_3arcminA.arf",
        "sxt-s_120210_ts02um_intallpxl.arf"]
fs_ids = ['pn-med', 'acisi_aimpt_cy17', 'aciss_aimpty_cy17', 'nustar_3arcminA',
    'sxt-s_120210_ts02um_intallpxl']


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_sloshing_return_events':
        params = [[a, r] for a, r in zip(arfs, rmfs)]
        metafunc.parametrize('arf, rmf', params, ids=fs_ids)

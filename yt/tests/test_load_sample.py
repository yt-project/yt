import logging
import os
import textwrap
from time import time

import numpy as np
import pytest

from yt.config import ytcfg
from yt.frontends.enzo.api import EnzoDataset
from yt.loaders import load_sample
from yt.sample_data.api import get_data_registry_table
from yt.testing import requires_module_pytest
from yt.utilities.logger import ytLogger


@pytest.fixture()
def tmp_data_dir(tmp_path):
    pre_test_data_dir = ytcfg["yt", "test_data_dir"]
    ytcfg.set("yt", "test_data_dir", str(tmp_path))

    yield tmp_path

    ytcfg.set("yt", "test_data_dir", pre_test_data_dir)


@pytest.fixture()
def capturable_logger(caplog):
    """
    This set the minimal conditions to make pytest's caplog fixture usable.
    """

    propagate = ytLogger.propagate
    ytLogger.propagate = True

    with caplog.at_level(logging.INFO):
        yield

    ytLogger.propagate = propagate


@requires_module_pytest("pandas", "pooch")
@pytest.mark.usefixtures("capturable_logger")
def test_load_sample_small_dataset(tmp_data_dir, caplog):
    start = time()
    ds = load_sample("ToroShockTube", progressbar=False, timeout=30)
    elapsed1 = time() - start
    assert isinstance(ds, EnzoDataset)

    text = textwrap.dedent(
        f"""
            'ToroShockTube' is not available locally. Looking up online.
            Downloading from https://yt-project.org/data/ToroShockTube.tar.gz
            Untaring downloaded file to '{str(tmp_data_dir)}'
            Parameters: current_time              = 0.2
            Parameters: domain_dimensions         = [100   1   1]
            Parameters: domain_left_edge          = [0. 0. 0.]
            Parameters: domain_right_edge         = [1. 1. 1.]
            Parameters: cosmological_simulation   = 0
        """
    ).strip("\n")
    expected = [("yt", 20, message) for message in text.split("\n")]
    assert caplog.record_tuples == expected

    caplog.clear()
    # loading a second time should not result in a download request
    start = time()
    ds2 = load_sample("ToroShockTube")
    elapsed2 = time() - start
    assert isinstance(ds2, EnzoDataset)

    # caching should make a significant runtime difference even in fast networks
    assert elapsed2 < elapsed1

    text = textwrap.dedent(
        f"""
            Sample dataset found in '{os.path.join(str(tmp_data_dir), 'ToroShockTube')}'
            Parameters: current_time              = 0.2
            Parameters: domain_dimensions         = [100   1   1]
            Parameters: domain_left_edge          = [0. 0. 0.]
            Parameters: domain_right_edge         = [1. 1. 1.]
            Parameters: cosmological_simulation   = 0
        """
    ).strip("\n")
    expected = [("yt", 20, message) for message in text.split("\n")]
    assert caplog.record_tuples == expected


@requires_module_pytest("pandas", "pooch")
@pytest.mark.usefixtures("capturable_logger")
def test_load_sample_timeout(tmp_data_dir, caplog):
    from requests.exceptions import ConnectTimeout, ReadTimeout

    # note that requests is a direct dependency to pooch,
    # so we don't need to mark it in a decorator.
    with pytest.raises((ConnectTimeout, ReadTimeout)):
        load_sample("IsolatedGalaxy", progressbar=False, timeout=0.00001)


@requires_module_pytest("pandas", "requests")
@pytest.mark.xfail(reason="Registry is currently incomplete.")
def test_registry_integrity():
    reg = get_data_registry_table()
    assert not any(reg.isna())


@pytest.fixture()
def sound_subreg():
    # this selection is needed because the full dataset is incomplete
    # so we test only the values that can be parsed
    reg = get_data_registry_table()
    return reg[["size", "byte_size"]][reg["size"].notna()]


@requires_module_pytest("pandas", "requests")
def test_registry_byte_size_dtype(sound_subreg):
    from pandas import Int64Dtype

    assert sound_subreg["byte_size"].dtype == Int64Dtype()


@requires_module_pytest("pandas", "requests")
def test_registry_byte_size_sign(sound_subreg):
    np.testing.assert_array_less(0, sound_subreg["byte_size"])


@requires_module_pytest("pandas", "requests")
def test_unknown_filename():
    fake_name = "these_are_not_the_files_your_looking_for"
    with pytest.raises(KeyError) as err:
        load_sample(fake_name)
        assert err.exc == f"Could not find '{fake_name}' in the registry."

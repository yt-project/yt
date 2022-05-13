import logging
import os
import re
import sys
import textwrap

import numpy as np
import pytest

from yt.config import ytcfg
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
@pytest.mark.parametrize(
    "fn, archive, exact_loc, class_",
    [
        (
            "ToroShockTube",
            "ToroShockTube.tar.gz",
            "ToroShockTube/DD0001/data0001",
            "EnzoDataset",
        ),
        (
            "KeplerianDisk",
            "KeplerianDisk.tar.gz",
            "KeplerianDisk/disk.out1.00000.athdf",
            "AthenaPPDataset",
        ),
        # trying with an exact name as well
        (
            "KeplerianDisk/disk.out1.00000.athdf",
            "KeplerianDisk.tar.gz",
            "KeplerianDisk/disk.out1.00000.athdf",
            "AthenaPPDataset",
        ),
        # check this special case because it relies on implementations
        # details in the AMRVAC frontend (using parfiles)
        # and could easily fail to load. See GH PR #3343
        (
            "rmi_dust_2d",
            "rmi_dust_2d.tar.gz",
            "rmi_dust_2d/output0001.dat",
            "AMRVACDataset",
        ),
    ],
)
def test_load_sample_small_dataset(
    fn, archive, exact_loc, class_: str, tmp_data_dir, caplog
):
    ds = load_sample(fn, progressbar=False, timeout=30)
    assert type(ds).__name__ == class_

    text = textwrap.dedent(
        f"""
            '{fn.replace('/', os.path.sep)}' is not available locally. Looking up online.
            Downloading from https://yt-project.org/data/{archive}
            Untaring downloaded file to '{str(tmp_data_dir)}'
        """
    ).strip("\n")
    expected = [("yt", 20, message) for message in text.split("\n")]
    assert caplog.record_tuples[:3] == expected

    caplog.clear()
    # loading a second time should not result in a download request
    ds2 = load_sample(fn)
    assert type(ds2).__name__ == class_

    assert caplog.record_tuples[0] == (
        "yt",
        20,
        f"Sample dataset found in '{os.path.join(str(tmp_data_dir), *exact_loc.split('/'))}'",
    )


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
    with pytest.raises(ValueError, match=f"'{fake_name}' is not an available dataset."):
        load_sample(fake_name)


@requires_module_pytest("pandas", "requests")
def test_typo_filename():
    wrong_name = "Isolatedgalaxy"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"'{wrong_name}' is not an available dataset. Did you mean 'IsolatedGalaxy' ?"
        ),
    ):
        load_sample(wrong_name, timeout=1)


@pytest.fixture()
def fake_data_dir_in_a_vaccum(tmp_path):
    pre_test_data_dir = ytcfg["yt", "test_data_dir"]
    ytcfg.set("yt", "test_data_dir", "/does/not/exist")
    origin = os.getcwd()
    os.chdir(tmp_path)

    yield

    ytcfg.set("yt", "test_data_dir", pre_test_data_dir)
    os.chdir(origin)


@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="can't figure out how to match the warning message in a cross-platform way",
)
@requires_module_pytest("pandas", "pooch")
@pytest.mark.usefixtures("fake_data_dir_in_a_vaccum")
def test_data_dir_broken():
    # check that load_sample still works if the test_data_dir
    # isn't properly set, in which case we should dl to the
    # cwd and issue a warning.
    msg = (
        r"Storage directory from yt config doesn't exist "
        r"\(currently set to '/does/not/exist'\)\. "
        r"Current working directory will be used instead\."
    )
    with pytest.warns(UserWarning, match=msg):
        load_sample("ToroShockTube")

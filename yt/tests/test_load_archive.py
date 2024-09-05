import os
import sys
import tarfile
import time

import pytest

from yt.config import ytcfg
from yt.loaders import load_archive
from yt.sample_data.api import _download_sample_data_file, get_data_registry_table
from yt.testing import requires_module_pytest
from yt.utilities.exceptions import YTUnidentifiedDataType


@pytest.fixture()
def data_registry():
    yield get_data_registry_table()


@pytest.fixture()
def tmp_data_dir(tmp_path):
    pre_test_data_dir = ytcfg["yt", "test_data_dir"]
    ytcfg.set("yt", "test_data_dir", str(tmp_path))

    yield tmp_path

    ytcfg.set("yt", "test_data_dir", pre_test_data_dir)


# Note: ratarmount cannot currently be installed on Windows as of v0.8.1
@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="ratarmount cannot currently be installed on Windows as of v0.8.1",
)
@pytest.mark.skipif(
    os.environ.get("JENKINS_HOME") is not None,
    reason="Archive mounting times out on Jenkins.",
)
@requires_module_pytest("pooch", "ratarmount")
@pytest.mark.parametrize(
    "fn, exact_loc, class_",
    [
        (
            "ToroShockTube.tar.gz",
            "ToroShockTube/DD0001/data0001",
            "EnzoDataset",
        ),
        (
            "ramses_sink_00016.tar.gz",
            "ramses_sink_00016/output_00016",
            "RAMSESDataset",
        ),
    ],
)
@pytest.mark.parametrize("archive_suffix", ["", ".gz"])
def test_load_archive(
    fn, exact_loc, class_: str, archive_suffix, tmp_data_dir, data_registry
):
    # Download the sample .tar.gz'd file
    targz_path = _download_sample_data_file(filename=fn)
    tar_path = targz_path.with_suffix(archive_suffix)

    if tar_path != targz_path:
        # Open the tarfile and uncompress it to .tar, .tar.gz, and .tar.bz2 files
        with tarfile.open(targz_path, mode="r:*") as targz:
            mode = "w" + archive_suffix.replace(".", ":")
            with tarfile.open(tar_path, mode=mode) as tar:
                for member in targz.getmembers():
                    content = targz.extractfile(member)
                    tar.addfile(member, fileobj=content)

    # Now try to open the .tar.* files
    warn_msg = "The 'load_archive' function is still experimental and may be unstable."
    with pytest.warns(UserWarning, match=warn_msg):
        ds = load_archive(tar_path, exact_loc, mount_timeout=10)
    assert type(ds).__name__ == class_

    # Make sure the index is readable
    ds.index

    # Check cleanup
    mount_path = tar_path.with_name(tar_path.name + ".mount")
    assert mount_path.is_mount()

    ## Manually dismount
    ds.dismount()

    ## The dismounting happens concurrently, wait a few sec.
    time.sleep(2)

    ## Mount path should not exist anymore *and* have been deleted
    assert not mount_path.is_mount()
    assert not mount_path.exists()


@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="ratarmount cannot currently be installed on Windows as of v0.8.1",
)
@pytest.mark.skipif(
    os.environ.get("JENKINS_HOME") is not None,
    reason="Archive mounting times out on Jenkins.",
)
@pytest.mark.filterwarnings(
    "ignore:The 'load_archive' function is still experimental and may be unstable."
)
@requires_module_pytest("pooch", "ratarmount")
def test_load_invalid_archive(tmp_data_dir, data_registry):
    # Archive does not exist
    with pytest.raises(FileNotFoundError):
        load_archive("this_file_does_not_exist.tar.gz", "invalid_location")

    targz_path = _download_sample_data_file(filename="ToroShockTube.tar.gz")
    # File does not exist
    with pytest.raises(FileNotFoundError):
        load_archive(targz_path, "invalid_location")

    # File exists but is not recognized
    with pytest.raises(YTUnidentifiedDataType):
        load_archive(targz_path, "ToroShockTube/DD0001/data0001.memorymap")

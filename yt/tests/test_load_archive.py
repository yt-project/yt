import tarfile
import time

import pytest

from yt.config import ytcfg
from yt.loaders import load_archive
from yt.sample_data.api import _download_sample_data_file, get_data_registry_table
from yt.testing import requires_module_pytest


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
@requires_module_pytest("pooch")
@requires_module_pytest("ratarmount")
@pytest.mark.parametrize(
    "fn ,exact_loc, class_",
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
def test_load_archive(fn, exact_loc, class_: str, tmp_data_dir, data_registry):
    # Download the sample .tar.gz'd file
    targz_path = _download_sample_data_file(filename=fn)
    tar_paths = [targz_path.with_suffix(suffix) for suffix in ("", ".xz", ".bz2")]

    # Open the tarfile and uncompress it to .tar, .tar.gz, and .tar.bz2 files
    with tarfile.open(targz_path, mode="r:*") as targz:
        for tar_path in tar_paths:
            mode = "w" + tar_path.suffix.replace(".", ":")
            with tarfile.open(tar_path, mode=mode) as tar:
                for member in targz.getmembers():
                    content = targz.extractfile(member)
                    tar.addfile(member, fileobj=content)

    # Now try to open the .tar.* files
    for path in [targz_path] + tar_paths:
        ds = load_archive(path, exact_loc, mount_timeout=10)
        assert type(ds).__name__ == class_

        # Make sure the index is readable
        ds.index

        # Check cleanup
        mount_path = path.with_name(path.name + ".mount")
        assert mount_path.is_mount()

        ## Manually dismount
        ds.dismount()

        ## The dismounting happens concurrently, wait a few sec.
        time.sleep(2)

        ## Mount path should not exist anymore *and* have been deleted
        assert not mount_path.is_mount()
        assert not mount_path.exists()

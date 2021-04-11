import os
import sys
from pathlib import Path

import pytest

from yt.config import ytcfg_defaults
from yt.testing import assert_raises

IDS = [key for key in ytcfg_defaults["yt"].keys() if key.endswith("_dir")]
DIR_PATHS = [ytcfg_defaults["yt"][k] for k in IDS]


@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="a different behaviour is expected on Windows",
)
@pytest.mark.parametrize("path", DIR_PATHS, ids=IDS)
def test_default_data_dirs_posix(path):
    p = Path(path)

    # guard assert: this is essential to the integrity of the following try block
    assert not p.exists()

    try:
        with assert_raises(OSError):
            # "[Errno 30] Read-only file system:"
            os.makedirs(p)
    finally:
        # this is okay because the try block is guarded by the previous assert
        if p.exists():
            os.removedirs(p)


@pytest.mark.skipif(
    not sys.platform.startswith("win"),
    reason="a different behaviour is expected on POSIX based OSes",
)
@pytest.mark.parametrize("path", DIR_PATHS, ids=IDS)
def test_default_data_dirs_windows(path):
    p = Path(path)
    with assert_raises(OSError):
        # "[WinError 123] The filename, directory name, or volume label syntax is incorrect:"
        p.exists()

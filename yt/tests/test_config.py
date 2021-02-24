import os
import platform
from pathlib import Path

from yt.config import ytcfg_defaults
from yt.testing import assert_raises


def test_default_data_dirs():
    for key, val in ytcfg_defaults["yt"].items():
        if not key.endswith("_dir"):
            continue

        p = Path(val)
        if platform.system() == "Windows":
            with assert_raises(OSError):
                # "[WinError 123] The filename, directory name, or volume label syntax is incorrect:"
                p.exists()
        else:
            # guard assert: this is essential to the integrity of the following try block
            assert not p.exists()

        try:
            with assert_raises(OSError):
                # "[Errno 30] Read-only file system:"
                os.makedirs(p)
        finally:
            # this is okay because the try block is guarded by
            # the previous assert
            if platform.system() != "Windows" and p.exists():
                os.removedirs(p)

import os
import shutil
import tempfile
from unittest import TestCase

import yt
from yt.testing import fake_random_ds


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


class SigmaClipTest(TestCase):
    # This toggles using a temporary directory. Turn off to examine images.
    use_tmpdir = True

    def setUp(self):
        if self.use_tmpdir:
            self.curdir = os.getcwd()
            # Perform I/O in safe place instead of yt main dir
            self.tmpdir = tempfile.mkdtemp()
            os.chdir(self.tmpdir)
        else:
            self.curdir, self.tmpdir = None, None

    def tearDown(self):
        if self.use_tmpdir:
            os.chdir(self.curdir)
            shutil.rmtree(self.tmpdir)

    def test_sigma_clip(self):
        ds = fake_random_ds(32)
        sc = yt.create_scene(ds)
        sc.save("clip_2.png", sigma_clip=2)

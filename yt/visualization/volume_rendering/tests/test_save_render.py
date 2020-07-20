import os
import tempfile
import shutil
import yt
from yt.testing import fake_random_ds
from unittest import TestCase

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"


class SaveRenderTest(TestCase):
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

    def test_save_render(self):
        ds = fake_random_ds(ndims=32)
        sc = yt.create_scene(ds)
        sc.save('raw.png') # will use render = True by default
        sc.save('clip_2.png', sigma_clip=2, render=False) # will pull render
        sc.save('clip_4.png', sigma_clip=4.0, render=False)

        return sc

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


class SaveRenderTest(TestCase):
    # This toggles using a temporary directory. Turn off to examine images.
    use_tmpdir = True
    tmpdir = "./"

    def setUp(self):
        if self.use_tmpdir:
            tempfile.mkdtemp()
            self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        if self.use_tmpdir:
            shutil.rmtree(self.tmpdir)

    def test_save_render(self):
        ds = fake_random_ds(ndims=32)
        sc = yt.create_scene(ds)

        # make sure it renders if nothing exists, even if render = False
        sc.save(os.path.join(self.tmpdir, "raw.png"), render=False)
        # make sure it re-renders
        sc.save(os.path.join(self.tmpdir, "raw_2.png"), render=True)
        # make sure sigma clip does not re-render
        sc.save(os.path.join(self.tmpdir, "clip_2.png"), sigma_clip=2.0, render=False)
        sc.save(os.path.join(self.tmpdir, "clip_4.png"), sigma_clip=4.0, render=False)

        # save a different format with/without sigma clips
        sc.save(os.path.join(self.tmpdir, "no_clip.jpg"), render=False)
        sc.save(os.path.join(self.tmpdir, "clip_2.jpg"), sigma_clip=2, render=False)

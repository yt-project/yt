"""
Tests for loading in-memory datasets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import shutil
import tempfile
import unittest

from yt.testing import assert_raises
from yt.convenience import load
from yt.utilities.exceptions import YTOutputNotIdentified

class TestEmptyLoad(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.curdir = os.getcwd()
        os.chdir(self.tmpdir)

        # create 0 byte file
        open("empty_file", "a")

        # create empty directory
        os.makedirs("empty_directory")

    def tearDown(self):
        os.chdir(self.curdir)
        shutil.rmtree(self.tmpdir)

    def test_load_empty_file(self):
        assert_raises(YTOutputNotIdentified, load, "not_a_file")
        assert_raises(YTOutputNotIdentified, load, "empty_file")
        assert_raises(YTOutputNotIdentified, load, "empty_directory")

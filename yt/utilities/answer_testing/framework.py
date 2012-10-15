"""
Answer Testing using Nose as a starting point

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import logging
import os
import shelve

from yt.testing import *
from yt.utilities.command_line import get_yt_version
from yt.config import ytcfg
from nose.plugins import Plugin

log = logging.getLogger('nose.plugins.answer-testing')

class AnswerTesting(Plugin):
    name = "answer-testing"

    def options(self, parser, env=os.environ):
        super(AnswerTesting, self).options(parser, env=env)
        test_storage_directory = ytcfg.get("yt", "test_storage_dir")
        try:
            my_hash = get_yt_version()
        except:
            my_hash = "UNKNOWN%s" % (time.time())
        parser.add_option("--answer-parameter-file", dest="parameter_file",
            default=os.path.join(os.getcwd(), "tests/DD0010/moving7_0010"),
            help="The parameter file value to feed to 'load' to test against")
        parser.add_option("--answer-output", dest="storage_dir",
            default=test_storage_directory,
            help="Base directory for storing test output.")
        parser.add_option("--answer-compare", dest="compare_name",
            default=None,
            help="The name against which we will compare")
        parser.add_option("--answer-name", dest="this_name",
            default=my_hash,
            help="The name we'll call this set of tests")

    def configure(self, options, conf):
        super(AnswerTesting, self).configure(options, conf)
        if not self.enabled:
            return
        AnswerTestingTest.result_storage = shelve.Shelf(
            os.path.join(options.storage_dir,
                         options.this_name))
        if options.compare_name is not None:
            AnswerTestingTest.reference_storage = shelve.Shelf(
                os.path.join(options.storage_dir,
                            options.compare_name))

    def finalize(self, result):
        pass

class AnswerTestingTest(object):
    reference_storage = None

    description = None
    def __init__(self, name, pf_fn):
        self.pf = load(pf_fn)
        self.name = "%s_%s" % (pf, name)

    def __call__(self):
        if self.reference_storage is not None:
            ov = self.reference_storage.get(self.name, None)
        else:
            ov = None
        nv = self.run()
        return self.compare(nv, ov)

    def compare(self, new_result, old_result):
        pass

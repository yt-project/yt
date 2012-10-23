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
import hashlib
import contextlib
import urllib2
import cPickle
from collections import defaultdict
from nose.plugins import Plugin

run_big_data = False

class AnswerTesting(Plugin):
    name = "answer-testing"

    def options(self, parser, env=os.environ):
        super(AnswerTesting, self).options(parser, env=env)
        parser.add_option("--answer-compare", dest="compare_name",
            default=None, help="The name against which we will compare")
        parser.add_option("--answer-big-data", dest="big_data",
            default=False, help="Should we run against big data, too?",
            action="store_true")
        parser.add_option("--answer-name", dest="this_name",
            default=None,
            help="The name we'll call this set of tests")
        parser.add_option("--answer-store", dest="store_results",
            default=False, action="store_true")

    def configure(self, options, conf):
        super(AnswerTesting, self).configure(options, conf)
        if not self.enabled:
            return
        from yt.utilities.logger import disable_stream_logging
        disable_stream_logging()
        from yt.utilities.command_line import get_yt_version
        try:
            my_hash = get_yt_version()
        except:
            my_hash = "UNKNOWN%s" % (time.time())
        if options.this_name is None: options.this_name = my_hash
        from yt.config import ytcfg
        ytcfg["yt","__withintesting"] = "True"
        from .framework import AnswerTestingTest, AnswerTestOpener
        AnswerTestingTest.result_storage = \
            self.result_storage = defaultdict(dict)
        if options.compare_name is not None:
            # Now we grab from our S3 store
            if options.compare_name == "latest":
                options.compare_name = _latest
            AnswerTestingTest.reference_storage = \
                AnswerTestOpener(options.compare_name)
        self.answer_name = options.this_name
        self.store_results = options.store_results
        global run_big_data
        run_big_data = options.big_data

    def finalize(self, result):
        # This is where we dump our result storage up to Amazon, if we are able
        # to.
        if self.store_results is False: return
        import boto
        from boto.s3.key import Key
        c = boto.connect_s3()
        bucket = c.get_bucket("yt-answer-tests")
        for pf_name in self.result_storage:
            rs = cPickle.dumps(self.result_storage[pf_name])
            tk = bucket.get_key("%s_%s" % (self.answer_name, pf_name)) 
            if tk is not None: tk.delete()
            k = Key(bucket)
            k.key = "%s_%s" % (self.answer_name, pf_name)
            k.set_contents_from_string(rs)
            k.set_acl("public-read")


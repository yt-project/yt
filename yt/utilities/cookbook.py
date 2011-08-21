"""
A way to find an utilize recipes

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

# See also:
#  http://en.wikipedia.org/wiki/Me_(mythology)

import os
import argparse
import abc
import glob
import imp
import types
import sys
import pprint

def _load_intent(intent_path):
    mname = os.path.basename(intent_path[:-3])
    f, filename, desc = imp.find_module(mname,
            [os.path.dirname(intent_path)])
    intent = imp.load_module(mname, f, filename, desc)
    for i in sorted(dir(intent)):
        obj = getattr(intent, i)
        if issubclass(obj, Intent) and \
           isinstance(obj.desc, types.StringTypes):
             return obj
    return None

def _find_cookbook_dir():
    yt_dest = os.environ.get("YT_DEST", None)
    if yt_dest is None:
        print "YT_DEST is not set!  Set it and try again."
        return False
    cookbook_dir = os.path.join(yt_dest,
        "src/yt-supplemental/yt-cookbook/intents")
    if not os.path.isdir(cookbook_dir):
        print "Cookbook does not contain 'intents' directory."
        print "Update with 'yt instinfo -u' and try again."
        print "(%s)" % cookbook_dir
        return False
    return cookbook_dir

class Intent(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, args):
        self.args = args
        if "help" in self.args:
            print
            print "The arguments to supply, in order:"
            print
            print self.help
            print
            sys.exit()

    @abc.abstractmethod
    def run(self):
        pass

    @abc.abstractproperty
    def desc(self): pass

    @abc.abstractproperty
    def help(self): pass

    @classmethod
    def list_intents(self):
        intents = []
        cookbook_dir = _find_cookbook_dir()
        if cookbook_dir is False: return 1
        for fn in sorted(glob.glob(os.path.join(cookbook_dir, "*"))):
            # We skim them, looking for the 'Intent' subclass
            if any(("(Intent):" in line for line in open(fn))):
                intents.append((os.path.basename(fn)[:-3],
                                _load_intent(fn)))
        print
        print "Found these Intents:"
        print "\n".join(("% 15s: %s" % (a, b.desc) for a, b in intents))
        print

    @classmethod
    def select_intent(self, intent_name):
        cookbook_dir = _find_cookbook_dir()
        intent = None
        for fn in glob.glob(os.path.join(cookbook_dir, "*")):
            if os.path.basename(fn)[:-3] == intent_name:
                intent = _load_intent(fn)
        return intent

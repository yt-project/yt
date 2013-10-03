"""
Very simple convenience function for importing all the modules, setting up
the namespace and getting the last argument on the command line.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import

#
# ALL IMPORTS GO HERE
#

# First module imports
import os
from yt import *

# This next item will handle most of the actual startup procedures, but it will
# also attempt to parse the command line and set up the global state of various
# operations.  The variable unparsed_args is not used internally but is
# provided as a convenience for users who wish to parse arguments in scripts.
# See http://lists.spacepope.org/pipermail/yt-dev-spacepope.org/2011-December/
#     001727.html
from .config import ytcfg
from .startup_tasks import unparsed_args

from yt.utilities.logger import level as __level
from yt.config import ytcfgDefaults
if __level >= int(ytcfgDefaults["loglevel"]):
    # This won't get displayed.
    mylog.debug("Turning off NumPy error reporting")
    np.seterr(all = 'ignore')

# We load plugins.  Keep in mind, this can be fairly dangerous -
# the primary purpose is to allow people to have a set of functions
# that get used every time that they don't have to *define* every time.
# This way, other command-line tools can be used very simply.
# Unfortunately, for now, I think the easiest and simplest way of doing
# this is also the most dangerous way.
if ytcfg.getboolean("yt","loadfieldplugins"):
    enable_plugins()

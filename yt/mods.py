#
# ALL IMPORTS GO HERE
#

import os

# This next item will handle most of the actual startup procedures, but it will
# also attempt to parse the command line and set up the global state of various
# operations.  The variable unparsed_args is not used internally but is
# provided as a convenience for users who wish to parse arguments in scripts.
# https://mail.python.org/archives/list/yt-dev@python.org/thread/L6AQPJ3OIMJC5SNKVM7CJG32YVQZRJWA/
import yt.startup_tasks as __startup_tasks
from yt import *
from yt.config import ytcfg, ytcfg_defaults
from yt.utilities.logger import level as __level

unparsed_args = __startup_tasks.unparsed_args


if __level >= int(ytcfg_defaults["loglevel"]):
    # This won't get displayed.
    mylog.debug("Turning off NumPy error reporting")
    np.seterr(all="ignore")

# We load plugins.  Keep in mind, this can be fairly dangerous -
# the primary purpose is to allow people to have a set of functions
# that get used every time that they don't have to *define* every time.
# This way, other command-line tools can be used very simply.
# Unfortunately, for now, I think the easiest and simplest way of doing
# this is also the most dangerous way.
if ytcfg.getboolean("yt", "loadfieldplugins"):
    enable_plugins()

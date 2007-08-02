"""
Fido is a module for moving and storing data.  It is designed to watch for
outputs, as well as transparent archiving them during and after a simulation
run.  Additionally, function handlers can be called to deal with output upon
its creation.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.logger import fidoLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

mylog.warning("shutil used in FileHandling may lead to sub-optimal performance")
import os, os.path, shutil, time, sys, glob, types

WAITBETWEEN=5
# This should go in ytcfg , but it isn't there just yet.
OTHERFILES=["rates.out","cool_rates.out"]
NEW_OUTPUT_CREATED = "newOutput"

# Let's define some useful functions.  Maybe these should go elsewhere?

def getParameterLine(filename, parameter):
    f = open(filename)
    lines = filter(lambda a: a.startswith(parameter), f)
    if len(lines) > 1:
        raise KeyError, "More than one line matches that parameter!"
    else: return lines[0]

def getParentDir(filename):
    return os.path.normpath( \
            os.path.split(os.path.dirname(os.path.abspath(filename)))[0])

from OutputCollection import *
from FileHandling import *
from OutputWatcher import *

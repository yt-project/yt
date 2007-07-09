"""


"""

from yt.logger import fidoLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

mylog.warning("shutil used in FileHandling may lead to sub-optimal performance")
import os, os.path, shutil, time, sys, glob

WAITBETWEEN=5
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

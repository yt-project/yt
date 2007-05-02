"""
Fido
====

    Fido is the messenger/protector of data.  It takes data outputs, puts them
    wherever you want, and then calls a function handler to deal with that data.
    Ultimately Fido will deal with all file-handling; submission of runs to a
    central (user-specific) database is in the works, and Fido will be the entity
    that moves the files in and out of storage.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

# Nothing here yet

from yt.logger import fidoLogger as mylog
from yt.config import ytcfg

import glob, os.path, time, ConfigParser
from filehandler import *
from dbhandler import *
from watch import *
from walkies import run, selectRun, selectPF, runBrowse, outputBrowse

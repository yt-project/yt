"""
Raven
=====

    Raven is the plotting interface, with support for several
    different engines.  Well, two for now, but maybe more later.
    Who knows?

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""
from yt.logger import ravenLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *
import yt.lagos as lagos
try:
    import yt.deliverator
except:
    mylog.warning("Deliverator import failed; all deliverator actions will fail!")

import time, types, string, os

from PlotConfig import *

# We now check with ytcfg to see which backend we want

backend = ytcfg["raven","backend"]

if backend.upper()=="HD":
    import backends.HD as be

from PlotTypes import *

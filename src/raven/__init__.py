"""
Raven
=====

    Raven is the  interface to HippoDraw.  All data plotting goes through
    Raven.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""
from yt.logger import ravenLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

import time, types, string, os

try:
    import hippo
except ImportError:
    mylog.warning("hippo import failed")
import yt.lagos as lagos

from EnzoPlotConfig import *
from EnzoPlotTypes import *
from EnzoHippoType import *
from EnzoGallery import *

try:
    import matplotlib
    from AMRPlot import *
except ImportError:
    mylog.warning("matplotlib failed to import")


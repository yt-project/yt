"""
Deliverator
===========

    The Deliverator is a TurboGears-based system for querying and displaying
    images.  Images are dumped from Raven into local, web-accessible storage space,
    and then metadata about those images is submitted to The Deliverator.  The user
    (you) then goes to the Deliverator website and views those plots.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import yt.lagos as lagos
from yt.config import ytcfg
from yt.logger import deliveratorLogger as mylog
#from services_server import *
from services_types import *
from services import *
from upload import *


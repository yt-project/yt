"""
Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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


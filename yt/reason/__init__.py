"""
Initializer for reason.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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


from yt.config import ytcfg
import yt.logger
mylog = yt.logger.reasonLogger

import yt.lagos as lagos
import yt.raven as raven
import yt.fido as fido

from yt.arraytypes import *

from math import log10, sqrt

import os, types, Toolbars

# Now let's set up the matplotlib stuff
import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
import pylab

#import yt.raven.backends.MPL as be
from yt.raven import be, color_maps
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
import wx, wx.py, wx.aui
import wx.lib.mixins.listctrl as listmix
from wx.lib.pubsub import Publisher
import matplotlib.backends.backend_wx as be_wx
import matplotlib.figure

from LoggingSetup import *
from Functions import *
from Windows import *
from Notebook import *
from App import *

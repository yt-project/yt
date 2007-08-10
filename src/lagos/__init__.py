"""
Lagos
=====

    Lagos defines a set of class structures for Enzo data handling.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

"""

from yt.logger import lagosLogger as mylog
from yt.config import ytcfg

try:
    from pyhdf_np import SD          # NumPy
    import pyhdf_np.error  # NumPy
except:
    mylog.warning("No HDF4 support")
import warnings, tables
warnings.simplefilter("ignore", tables.NaturalNameWarning)

from yt.arraytypes import *
import weakref
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os, os.path, types, exceptions, re
from stat import ST_CTIME

#import RavenCombine, fields, chemistry
import time

try:
    from yt.enki import EnzoInterface
except:
    pass

try:
    import EnzoFortranRoutines
except:
    pass
    #mylog.warning("EnzoFortranRoutines import failed; all fortran calls will fail!")

import PointCombine
from EnzoDefs import *
from DerivedFields import fieldInfo, log_fields, colormap_dict
from DerivedQuantities import quantityInfo
from EnzoFortranWrapper import cosmologyGetUnits
from DataReadingFuncs import *
from BaseGridType import *
from ClusterFiles import *
from BaseDataTypes import *
from EnzoRateData import *
from HierarchyType import *
from OutputTypes import *
from EnzoRunType import *
#from ThreadsAndQueues import *
from DataCube import *

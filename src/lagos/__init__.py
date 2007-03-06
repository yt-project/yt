"""
Lagos
=====

    Lagos defines a set of class structures for Enzo data handling.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

"""


from pyhdf import SD
import pyhdf.error
import tables, warnings
from numarray import *
import numarray.objects as obj
import numarray.nd_image as nd
import numarray.records as rec
import numarray
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os, os.path, types, exceptions, re
from stat import ST_CTIME

#import RavenCombine, fields, chemistry
import time

try:
    import EnzoFortranRoutines
except:
    mylog.warning("EnzoFortranRoutines import failed; all fortran calls will fail!")
from EnzoDefs import *
from EnzoDerivedFields import *
from EnzoFortranWrapper import *
from EnzoDataFuncs import *
from EnzoGridType import *
from EnzoDataTypes import *
from EnzoHierarchyType import *
from EnzoRunType import *
import EnzoCombine
from EnzoRateData import *

try:
    from yt.enki import EnzoInterface
except:
    pass

from yt.logger import lagosLogger as mylog


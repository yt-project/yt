"""
Enki
====

Enki is the package used to create data, and instantiate runs.  It supports
creating Enzo Problems, and then using SWIG-interfaced Enzo calls to actually
create the data for those problems.  Additionally, facilities are being
developed to use Enki to directly execute runs.

G{packagetree}

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import sys
from yt.logger import enkiLogger as mylog
from yt.config import ytcfg
from yt.arraytypes import *

# Now we import the SWIG enzo interface
# Note that we're going to try super-hard to get the one that's local for the
# user

sp = sys.path

if ytcfg.has_option("SWIG", "EnzoInterfacePath"):
    swig_path = ytcfg.get("SWIG","EnzoInterfacePath")
    mylog.info("Using %s as path to SWIG Interface", swig_path)
    sys.path = sys.path[:1] + [swig_path] + sys.path[1:] # We want '' to be the first

try:
    import EnzoInterface
    mylog.debug("Imported EnzoInterface successfully")
    has_SWIG = True
except ImportError, e:
    mylog.warning("EnzoInterface failed to import; all SWIG actions will fail")
    mylog.warning("(%s)", e)
    has_SWIG = False

sys.path = sp

from EnkiDefs import *

from HostProfiles import *
from CreateNewProblem import *
import mes
from Recompile import *

#from yt.lagos import * # Do we need this?
import yt.lagos as lagos

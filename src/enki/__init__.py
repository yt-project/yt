import sys
from yt.logger import enkiLogger as mylog
from yt.config import ytcfg

from numarray import *

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
except ImportError, e:
    mylog.warning("EnzoInterface failed to import; all SWIG actions will fail")
    mylog.warning("(%s)", e)

sys.path = sp

from HostProfiles import *
from CreateNewProblem import *
from EnzoProblem import *
import mes

from yt.lagos import *

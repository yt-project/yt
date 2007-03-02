from HostProfiles import *
import sys

from yt.logger import enkiLogger as mylog
from yt.config import ytcfg

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
except:
    mylog.warning("EnzoInterface failed to import; all SWIG actions will fail")

sys.path = sp

from yt.lagos import *

"""
Logging facility for yt
Will initialize everything, and associate one with each module

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import logging
import logging.handlers as handlers

logging.basicConfig(
    format="%(name)-10s %(levelname)-10s %(asctime)s %(message)s",
    level=logging.DEBUG
)

f = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s")

# We initialize each module individually
# Note that there are LOTS of ways we can log things, and some of these ought
# to be tossed into the home directory.
#
# Note that you can override the default logging in scripts that import yt,
# because all calls to logging.getLogger((WHATEVER) return the same logger
# object for a given value of WHATEVER.
#
# Note that if you don't have permissions to the current directory, it just
# dumps to /dev/null .  This may be undesirable.

try:
    fidoHandler = handlers.RotatingFileHandler("fido.log", maxBytes=10*1024*1024, backupCount=10)
except IOError:
    fidoHandler = logging.FileHandler("/dev/null")
fidoHandler.setFormatter(f)
fidoLogger = logging.getLogger("yt.fido")
fidoLogger.addHandler(fidoHandler)

try:
    ravenHandler = handlers.RotatingFileHandler("raven.log", maxBytes=10*1024*1024, backupCount=10)
except IOError:
    ravenHandler = logging.FileHandler("/dev/null")
ravenHandler.setFormatter(f)
ravenLogger = logging.getLogger("yt.raven")
ravenLogger.addHandler(ravenHandler)
ravenLogger.setLevel(logging.NOTSET)

try:
    lagosHandler = handlers.RotatingFileHandler("lagos.log", maxBytes=10*1024*1024, backupCount=10)
except IOError:
    lagosHandler = logging.FileHandler("/dev/null")
lagosHandler.setFormatter(f)
lagosLogger = logging.getLogger("yt.lagos")
lagosLogger.addHandler(lagosHandler)
lagosLogger.setLevel(logging.NOTSET)

try:
    enkiHandler = handlers.RotatingFileHandler("enki.log", maxBytes=10*1024*1024, backupCount=10)
except IOError:
    enkiHandler = logging.FileHandler("/dev/null")
enkiHandler.setFormatter(f)
enkiLogger = logging.getLogger("yt.enki")
enkiLogger.addHandler(enkiHandler)
enkiLogger.setLevel(logging.NOTSET)

try:
    deliveratorHandler = handlers.RotatingFileHandler("deliverator.log", maxBytes=10*1024*1024, backupCount=10)
except IOError:
    deliveratorHandler = logging.FileHandler("/dev/null")
deliveratorHandler.setFormatter(f)
deliveratorLogger = logging.getLogger("yt.deliverator")
deliveratorLogger.addHandler(deliveratorHandler)
deliveratorLogger.setLevel(logging.NOTSET)

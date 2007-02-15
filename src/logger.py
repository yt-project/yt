# Logging facility for yt
# Will initialize everything, and associate one with each module

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

fidoHandler = handlers.RotatingFileHandler("fido.log", maxBytes=10*1024, backupCount=10)
fidoHandler.setFormatter(f)
fidoLogger = logging.getLogger("yt.fido")
fidoLogger.addHandler(fidoHandler)

ravenHandler = handlers.RotatingFileHandler("raven.log", maxBytes=10*1024, backupCount=10)
ravenHandler.setFormatter(f)
ravenLogger = logging.getLogger("yt.raven")
ravenLogger.addHandler(ravenHandler)
ravenLogger.setLevel(logging.NOTSET)

lagosHandler = handlers.RotatingFileHandler("lagos.log", maxBytes=10*1024, backupCount=10)
lagosHandler.setFormatter(f)
lagosLogger = logging.getLogger("yt.lagos")
lagosLogger.addHandler(lagosHandler)
lagosLogger.setLevel(logging.NOTSET)

deliveratorHandler = handlers.RotatingFileHandler("deliverator.log", maxBytes=10*1024, backupCount=10)
deliveratorHandler.setFormatter(f)
deliveratorLogger = logging.getLogger("yt.deliverator")
deliveratorLogger.addHandler(deliveratorHandler)
deliveratorLogger.setLevel(logging.NOTSET)

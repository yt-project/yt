"""
Logging facility for yt
Will initialize everything, and associate one with each module

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import logging, os
import logging.handlers as handlers
from yt.config import ytcfg

level = min(max(ytcfg.getint("yt","loglevel"),0),50)
fstring = "%(name)-10s %(levelname)-10s %(asctime)s %(message)s"
logging.basicConfig(
    format=fstring,
    level=level
)

f = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s")

ytLogger = logging.getLogger("yt")
ytLogger.debug("Set log level to %s", level)

fidoLogger = logging.getLogger("yt.fido")
ravenLogger = logging.getLogger("yt.raven")
lagosLogger = logging.getLogger("yt.lagos")
enkiLogger = logging.getLogger("yt.enki")
deliveratorLogger = logging.getLogger("yt.deliverator")

# Maybe some day we'll make this more configurable...  unfortunately, for now,
# we preserve thread-safety by opening in the current directory.

mb = 10*1024*1024
bc = 10

if ytcfg.getboolean("yt","logfile") and os.access(".", os.W_OK):
    if ytcfg.getboolean("yt","unifiedlogfile"):
        ytHandler = handlers.RotatingFileHandler("yt.log", maxBytes=mb, backupCount=bc)
        k = logging.Formatter(fstring)
        ytHandler.setFormatter(k)
        ytLogger.addHandler(ytHandler)
        fidoLogger.addHandler(ytHandler)
        ravenLogger.addHandler(ytHandler)
        lagosLogger.addHandler(ytHandler)
        enkiLogger.addHandler(ytHandler)
        deliveratorLogger.addHandler(ytHandler)
    else:  # If we *don't* want a unified file handler (which is the default now!)
        fidoHandler = handlers.RotatingFileHandler("fido.log", maxBytes=mb, backupCount=bc)
        fidoHandler.setFormatter(f)
        fidoLogger.addHandler(fidoHandler)

        ravenHandler = handlers.RotatingFileHandler("raven.log", maxBytes=mb, backupCount=bc)
        ravenHandler.setFormatter(f)
        ravenLogger.addHandler(ravenHandler)

        lagosHandler = handlers.RotatingFileHandler("lagos.log", maxBytes=mb, backupCount=bc)
        lagosHandler.setFormatter(f)
        lagosLogger.addHandler(lagosHandler)

        enkiHandler = handlers.RotatingFileHandler("enki.log", maxBytes=mb, backupCount=bc)
        enkiHandler.setFormatter(f)
        enkiLogger.addHandler(enkiHandler)

        deliveratorHandler = handlers.RotatingFileHandler("deliverator.log", maxBytes=mb, backupCount=bc)
        deliveratorHandler.setFormatter(f)
        deliveratorLogger.addHandler(deliveratorHandler)

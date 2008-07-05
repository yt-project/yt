"""
Logging facility for yt
Will initialize everything, and associate one with each module

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

import logging, os
import logging.handlers as handlers
from yt.config import ytcfg

level = min(max(ytcfg.getint("yt", "loglevel"), 0), 50)
fstring = "%(name)-10s %(levelname)-10s %(asctime)s %(message)s"
logging.basicConfig(
    format=fstring,
    level=level
)

f = logging.Formatter("%(levelname)-10s %(asctime)s %(message)s")

rootLogger = logging.getLogger()

ytLogger = logging.getLogger("yt")
ytLogger.debug("Set log level to %s", level)

fidoLogger = logging.getLogger("yt.fido")
ravenLogger = logging.getLogger("yt.raven")
lagosLogger = logging.getLogger("yt.lagos")
enkiLogger = logging.getLogger("yt.enki")
deliveratorLogger = logging.getLogger("yt.deliverator")
reasonLogger = logging.getLogger("yt.reason")

# Maybe some day we'll make this more configurable...  unfortunately, for now,
# we preserve thread-safety by opening in the current directory.

mb = 10*1024*1024
bc = 10

if ytcfg.getboolean("yt","logfile") and os.access(".", os.W_OK):
    if ytcfg.getboolean("yt","unifiedlogfile"):
        log_file_name = ytcfg.get("yt","LogFileName")
        ytHandler = handlers.RotatingFileHandler(log_file_name,
                                                 maxBytes=mb, backupCount=bc)
        k = logging.Formatter(fstring)
        ytHandler.setFormatter(k)
        ytLogger.addHandler(ytHandler)
    else:
        # If we *don't* want a unified file handler (which is the default now!)
        fidoHandler = handlers.RotatingFileHandler("fido.log",
                                                   maxBytes=mb, backupCount=bc)
        fidoHandler.setFormatter(f)
        fidoLogger.addHandler(fidoHandler)

        ravenHandler = handlers.RotatingFileHandler("raven.log",
                                                    maxBytes=mb, backupCount=bc)
        ravenHandler.setFormatter(f)
        ravenLogger.addHandler(ravenHandler)

        lagosHandler = handlers.RotatingFileHandler("lagos.log",
                                                    maxBytes=mb, backupCount=bc)
        lagosHandler.setFormatter(f)
        lagosLogger.addHandler(lagosHandler)

        enkiHandler = handlers.RotatingFileHandler("enki.log",
                                                   maxBytes=mb, backupCount=bc)
        enkiHandler.setFormatter(f)
        enkiLogger.addHandler(enkiHandler)

        deliveratorHandler = handlers.RotatingFileHandler("deliverator.log",
                                                maxBytes=mb, backupCount=bc)
        deliveratorHandler.setFormatter(f)
        deliveratorLogger.addHandler(deliveratorHandler)

        reasonHandler = handlers.RotatingFileHandler("reason.log",
                                                    maxBytes=mb, backupCount=bc)
        reasonHandler.setFormatter(f)
        reasonLogger.addHandler(reasonHandler)

def disable_stream_logging():
    # We just remove the root logger's handlers
    for handler in rootLogger.handlers:
        if isinstance(handler, logging.StreamHandler):
            rootLogger.removeHandler(handler)

if ytcfg.getboolean("yt","suppressStreamLogging"):
    disable_stream_logging()

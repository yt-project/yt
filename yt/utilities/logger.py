"""
Logging facility for yt
Will initialize everything, and associate one with each module

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2013 Matthew Turk.  All Rights Reserved.

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

import logging
from yt.config import ytcfg

# This next bit is grabbed from:
# http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored


def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if(levelno >= 50):
            color = '\x1b[31m'  # red
        elif(levelno >= 40):
            color = '\x1b[31m'  # red
        elif(levelno >= 30):
            color = '\x1b[33m'  # yellow
        elif(levelno >= 20):
            color = '\x1b[32m'  # green
        elif(levelno >= 10):
            color = '\x1b[35m'  # pink
        else:
            color = '\x1b[0m'  # normal
        ln = color + args[1].levelname + '\x1b[0m'
        args[1].levelname = ln
        return fn(*args)
    return new

level = min(max(ytcfg.getint("yt", "loglevel"), 0), 50)
ufstring = "%(name)-3s: [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s: [%(levelname)-18s] %(asctime)s %(message)s"
logging.basicConfig(
    format=ufstring,
    level=level
)

rootLogger = logging.getLogger()
ytLogger = logging.getLogger("yt")


def disable_stream_logging():
    # We just remove the root logger's handlers
    for handler in rootLogger.handlers:
        if isinstance(handler, logging.StreamHandler):
            rootLogger.removeHandler(handler)
    h = logging.NullHandler()
    ytLogger.addHandler(h)

original_emitter = logging.StreamHandler.emit


def colorize_logging():
    f = logging.Formatter(cfstring)
    if len(rootLogger.handlers) > 0:
        rootLogger.handlers[0].setFormatter(f)
    logging.StreamHandler.emit = add_coloring_to_emit_ansi(
        logging.StreamHandler.emit)


def uncolorize_logging():
    f = logging.Formatter(ufstring)
    if len(rootLogger.handlers) > 0:
        rootLogger.handlers[0].setFormatter(f)
    logging.StreamHandler.emit = original_emitter

if ytcfg.getboolean("yt", "coloredlogs"):
    colorize_logging()

if ytcfg.getboolean("yt", "suppressStreamLogging"):
    disable_stream_logging()

ytLogger.debug("Set log level to %s", level)

"""
Logging facility for yt
Will initialize everything, and associate one with each module



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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

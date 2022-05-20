import logging
import sys
from typing import Callable, Optional

from yt.utilities.configure import YTConfig, configuration_callbacks

_yt_sh: Optional[logging.StreamHandler] = None
_original_emitter: Optional[Callable[[logging.LogRecord], None]] = None


def set_log_level(level):
    """
    Select which minimal logging level should be displayed.

    Parameters
    ----------
    level: int or str
        Possible values by increasing level:
        0 or "notset"
        1 or "all"
        10 or "debug"
        20 or "info"
        30 or "warning"
        40 or "error"
        50 or "critical"
    """
    # this is a user-facing interface to avoid importing from yt.utilities in user code.

    if isinstance(level, str):
        level = level.upper()

    if level == "ALL":  # non-standard alias
        level = 1
    ytLogger.setLevel(level)
    ytLogger.debug("Set log level to %s", level)


ytLogger = logging.getLogger("yt")


class DuplicateFilter(logging.Filter):
    """A filter that removes duplicated successive log entries."""

    # source
    # https://stackoverflow.com/questions/44691558/suppress-multiple-messages-with-same-content-in-python-logging-module-aka-log-co
    def filter(self, record):
        current_log = (record.module, record.levelno, record.msg, record.args)
        if current_log != getattr(self, "last_log", None):
            self.last_log = current_log
            return True
        return False


ytLogger.addFilter(DuplicateFilter())


class DeprecatedFieldFilter(logging.Filter):
    """A filter that suppresses repeated logging of deprecated field warnings"""

    def __init__(self, name=""):
        self.logged_fields = []
        super().__init__(name=name)

    def filter(self, record):
        if not record.msg.startswith("The Derived Field"):
            return True

        field = record.args[0]
        if field in self.logged_fields:
            return False

        self.logged_fields.append(field)
        return True


ytLogger.addFilter(DeprecatedFieldFilter())

# This next bit is grabbed from:
# http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored


def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[0].levelno
        if levelno >= 50:
            color = "\x1b[31m"  # red
        elif levelno >= 40:
            color = "\x1b[31m"  # red
        elif levelno >= 30:
            color = "\x1b[33m"  # yellow
        elif levelno >= 20:
            color = "\x1b[32m"  # green
        elif levelno >= 10:
            color = "\x1b[35m"  # pink
        else:
            color = "\x1b[0m"  # normal
        ln = color + args[0].levelname + "\x1b[0m"
        args[0].levelname = ln
        return fn(*args)

    return new


ufstring = "%(name)-3s: [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s: [%(levelname)-18s] %(asctime)s %(message)s"


def colorize_logging():
    f = logging.Formatter(cfstring)
    ytLogger.handlers[0].setFormatter(f)
    ytLogger.handlers[0].emit = add_coloring_to_emit_ansi(ytLogger.handlers[0].emit)


def uncolorize_logging():
    global _original_emitter, _yt_sh
    if None not in (_original_emitter, _yt_sh):
        f = logging.Formatter(ufstring)
        ytLogger.handlers[0].setFormatter(f)
        _yt_sh.emit = _original_emitter


def disable_stream_logging():
    if len(ytLogger.handlers) > 0:
        ytLogger.removeHandler(ytLogger.handlers[0])
    h = logging.NullHandler()
    ytLogger.addHandler(h)


def _runtime_configuration(ytcfg: YTConfig) -> None:
    # only run this at the end of yt.__init__, after yt.config.ytcfg was instanciated

    global _original_emitter, _yt_sh

    if ytcfg.get("yt", "stdout_stream_logging"):
        stream = sys.stdout
    else:
        stream = sys.stderr

    _level = min(max(ytcfg.get("yt", "log_level"), 0), 50)

    if ytcfg.get("yt", "suppress_stream_logging"):
        disable_stream_logging()
    else:
        _yt_sh = logging.StreamHandler(stream=stream)
        # create formatter and add it to the handlers
        formatter = logging.Formatter(ufstring)
        _yt_sh.setFormatter(formatter)
        # add the handler to the logger
        ytLogger.addHandler(_yt_sh)
        ytLogger.setLevel(_level)
        ytLogger.propagate = False

        _original_emitter = _yt_sh.emit

        if ytcfg.get("yt", "colored_logs"):
            colorize_logging()


configuration_callbacks.append(_runtime_configuration)

import logging
import sys
from typing import Union

import pkg_resources
from rich.console import Console
from rich.logging import RichHandler
from rich.theme import Theme

from yt.config import loglevel_int2str, loglevel_str2int, ytcfg
from yt.utilities.exceptions import YTConfigError

# YTEP-0039: legacy helper functions to (de)activate color logging
# and suppress logging. They should be removed at some point.
# (see the YTEP)

ufstring = "%(name)-3s: [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s: [%(levelname)-18s] %(asctime)s %(message)s"

# This next bit is grabbed from:
# http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored


def add_coloring_to_emit_ansi(fn):
    # YTEP-0039: legacy only
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


def disable_stream_logging():
    # YTEP-0039: legacy only
    if len(ytLogger.handlers) > 0:
        ytLogger.removeHandler(ytLogger.handlers[0])
    h = logging.NullHandler()
    ytLogger.addHandler(h)


def colorize_logging():
    # YTEP-0039: legacy only
    f = logging.Formatter(cfstring)
    ytLogger.handlers[0].setFormatter(f)
    handler.emit = add_coloring_to_emit_ansi(handler.emit)


def uncolorize_logging():
    # YTEP-0039: legacy only
    try:
        f = logging.Formatter(ufstring)
        ytLogger.handlers[0].setFormatter(f)
        handler.emit = original_emitter
    except NameError:
        # `handler` and `original_emitter` are not defined because
        # ytcfg["logging", "stream"] == "none"
        # so we continue since there is nothing to uncolorize.
        pass


### end legacy-only definitions (YTEP-0039)


def set_log_level(level: Union[int, str]) -> None:
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

    if isinstance(level, int):
        ytcfg["logging", "level"] = loglevel_int2str(level)
    elif isinstance(level, str):
        ytcfg["logging", "level"] = level.upper()
        level = loglevel_str2int(level)
    else:
        raise TypeError(
            f"Expected an int or an str, got `{level}` with type `{type(level)}`."
        )

    ytLogger.setLevel(level)
    ytLogger.debug("log level set to %s", level)


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


ytLogger = logging.getLogger("yt")
ytLogger.addFilter(DuplicateFilter())

# normalize logging level string to all caps
ytcfg["logging", "level"] = ytcfg.get("logging", "level").upper()


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


_stream = ytcfg.get("logging", "stream").lower()
if _stream == "stdout":
    stream = sys.stdout
elif _stream == "stderr":
    stream = sys.stderr
elif _stream == "none":
    disable_stream_logging()
else:
    raise YTConfigError("yt.logging.stream", choices=["stdout", "stderr", "none"])

if _stream != "none":
    _handler = ytcfg.get("logging", "handler")

    if _handler == "legacy":
        handler = logging.StreamHandler(stream=stream)
        _fmt = ufstring
        _datefmt = None  # unspecified
    elif _handler == "rich":
        _width = ytcfg.get("logging", "width")
        width = None if _width <= 0 else _width

        theme_file = None
        no_color = not ytcfg.get("logging", "use_color")
        if no_color:
            # note that the monochrome theme doesn't deactivate rich markup (bold, italic, dim)
            theme_file = pkg_resources.resource_filename(
                "yt", "utilities/monochrome_logger_theme.ini"
            )

        _custom_theme = ytcfg.get("logging", "custom_theme")
        if _custom_theme:
            theme_file = _custom_theme

        if theme_file is None:
            theme = None
        else:
            theme = Theme.read(theme_file)

        console = Console(file=stream, width=width, theme=theme, no_color=no_color)
        handler = RichHandler(console=console)
        _fmt = ytcfg.get("logging", "format")
        _datefmt = ytcfg.get("logging", "date_format")
    else:
        raise YTConfigError("yt.logging.handler", choices=["legacy", "rich"])

    handler.setFormatter(logging.Formatter(fmt=_fmt, datefmt=_datefmt))

    # add the handler to the logger
    ytLogger.addHandler(handler)

    set_log_level(ytcfg.get("logging", "level"))
    ytLogger.propagate = False

    original_emitter = handler.emit
    if _handler == "legacy" and ytcfg.get("logging", "use_color"):
        colorize_logging()

# this exists only to support yt.mods
_level = loglevel_str2int(ytcfg.get("logging", "level"))

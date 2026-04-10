import io
import logging
import sys

from yt.utilities import logger as logger_mod


def _log_record(msg, levelno=20, args=(), module="test_module"):
    return logging.LogRecord(
        name="yt",
        level=levelno,
        pathname=__file__,
        lineno=1,
        msg=msg,
        args=args,
        exc_info=None,
    )


def test_duplicate_filter_and_deprecated_field_filter_behaviors():
    duplicate_filter = logger_mod.DuplicateFilter()
    first = _log_record("repeat", levelno=30, args=("x",))
    second = _log_record("repeat", levelno=30, args=("x",))
    third = _log_record("repeat", levelno=30, args=("y",))

    assert duplicate_filter.filter(first) is True
    assert duplicate_filter.filter(second) is False
    assert duplicate_filter.filter(third) is True

    deprecated_filter = logger_mod.DeprecatedFieldFilter()
    neutral = _log_record("ordinary message", args=("density",))
    deprecated = _log_record("The Derived Field %s is deprecated", args=("density",))

    assert deprecated_filter.filter(neutral) is True
    assert deprecated_filter.filter(deprecated) is True
    assert deprecated_filter.filter(deprecated) is False


def test_coloring_helpers_and_runtime_configuration_paths():
    original_handlers = list(logger_mod.ytLogger.handlers)
    original_level = logger_mod.ytLogger.level
    original_propagate = logger_mod.ytLogger.propagate
    original_emitter = logger_mod._original_emitter
    original_stream_handler = logger_mod._yt_sh

    try:
        handler = logging.StreamHandler(stream=io.StringIO())
        logger_mod.ytLogger.handlers = [handler]
        logger_mod._yt_sh = handler
        logger_mod._original_emitter = handler.emit

        logger_mod.colorize_logging()
        assert logger_mod.ytLogger.handlers[0].formatter._fmt == logger_mod.cfstring
        wrapped_emit = logger_mod.ytLogger.handlers[0].emit
        assert wrapped_emit is not logger_mod._original_emitter

        color_cases = {
            60: "\x1b[31m",
            40: "\x1b[31m",
            30: "\x1b[33m",
            20: "\x1b[32m",
            10: "\x1b[35m",
            0: "\x1b[0m",
        }
        for levelno, prefix in color_cases.items():
            record = _log_record("message", levelno=levelno)
            record.levelname = f"L{levelno}"
            wrapped_emit(record)
            assert record.levelname.startswith(prefix)

        logger_mod.uncolorize_logging()
        assert logger_mod.ytLogger.handlers[0].formatter._fmt == logger_mod.ufstring
        assert logger_mod._yt_sh.emit == logger_mod._original_emitter

        logger_mod._original_emitter = None
        logger_mod.uncolorize_logging()
        assert logger_mod.ytLogger.handlers[0].formatter._fmt == logger_mod.ufstring

        logger_mod.disable_stream_logging()
        assert isinstance(logger_mod.ytLogger.handlers[0], logging.NullHandler)

        class _FakeConfig:
            def __init__(self, values):
                self.values = values

            def get(self, section, key):
                return self.values[(section, key)]

        logger_mod.ytLogger.handlers = []
        logger_mod._original_emitter = None
        logger_mod._yt_sh = None
        logger_mod._runtime_configuration(
            _FakeConfig(
                {
                    ("yt", "stdout_stream_logging"): False,
                    ("yt", "log_level"): 12,
                    ("yt", "suppress_stream_logging"): True,
                    ("yt", "colored_logs"): False,
                }
            )
        )
        assert isinstance(logger_mod.ytLogger.handlers[0], logging.NullHandler)

        logger_mod.ytLogger.handlers = []
        logger_mod._original_emitter = None
        logger_mod._yt_sh = None
        logger_mod._runtime_configuration(
            _FakeConfig(
                {
                    ("yt", "stdout_stream_logging"): True,
                    ("yt", "log_level"): 99,
                    ("yt", "suppress_stream_logging"): False,
                    ("yt", "colored_logs"): True,
                }
            )
        )
        configured_handler = logger_mod.ytLogger.handlers[0]
        assert configured_handler.stream is sys.stdout
        assert logger_mod.ytLogger.level == 50
        assert logger_mod.ytLogger.propagate is False
        assert logger_mod._original_emitter is not None
        assert configured_handler.formatter._fmt == logger_mod.cfstring
    finally:
        logger_mod.ytLogger.handlers = original_handlers
        logger_mod.ytLogger.level = original_level
        logger_mod.ytLogger.propagate = original_propagate
        logger_mod._original_emitter = original_emitter
        logger_mod._yt_sh = original_stream_handler

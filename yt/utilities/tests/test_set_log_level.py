from yt.testing import assert_raises
from yt.utilities.logger import set_log_level


def test_valid_level():
    # test a subset of valid entries to cover
    # - case-insensitivity
    # - integer values
    # - "all" alias, which isn't standard
    for lvl in ("all", "ALL", 10, 42, "info", "warning", "ERROR", "CRITICAL"):
        set_log_level(lvl)


def test_invalid_level():
    # these are the exceptions raised by logging.Logger.setLog
    # since they are perfectly clear and readable, we check that nothing else
    # happens in the wrapper
    assert_raises(TypeError, set_log_level, 1.5)
    assert_raises(ValueError, set_log_level, "invalid_level")

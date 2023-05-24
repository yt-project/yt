# A namespace to store simple backports from most recent Python versions
# Backports should be contained in (if/else) blocks, checking on the runtime
# Python version, e.g.
#
# if sys.version_info < (3, 8):
#     ... # insert backported definitions
# else:
#     pass
#
# while the else part if meant as a no-op, it is required to allow automated
# cleanups with pyupgrade, when our minimal supported Python version reaches
# the requirement for backports to be unneeded.
#
# Likewise, imports from this module should be guarded by a similar condition,
# e.g.,
#
# if sys.version_info >= (3, 8):
#     from functools import cached_property
# else:
#     from yt._maintenance.backports import cached_property
import sys

if sys.version_info >= (3, 11):
    pass
else:
    from enum import Enum

    # backported from Python 3.11.0
    class ReprEnum(Enum):
        """
        Only changes the repr(), leaving str() and format() to the mixed-in type.
        """

    # backported from Python 3.11.0
    class StrEnum(str, ReprEnum):
        """
        Enum where members are also (and must be) strings
        """

        def __new__(cls, *values):
            "values must already be of type `str`"
            if len(values) > 3:
                raise TypeError(f"too many arguments for str(): {values!r}")
            if len(values) == 1:
                # it must be a string
                if not isinstance(values[0], str):
                    raise TypeError(f"{values[0]!r} is not a string")
            if len(values) >= 2:
                # check that encoding argument is a string
                if not isinstance(values[1], str):
                    raise TypeError(f"encoding must be a string, not {values[1]!r}")
            if len(values) == 3:
                # check that errors argument is a string
                if not isinstance(values[2], str):
                    raise TypeError("errors must be a string, not %r" % (values[2]))
            value = str(*values)
            member = str.__new__(cls, value)
            member._value_ = value
            return member

        def _generate_next_value_(name, start, count, last_values):  # noqa B902
            """
            Return the lower-cased version of the member name.
            """
            return name.lower()


builtin_zip = zip
if sys.version_info >= (3, 10):
    zip = builtin_zip
else:
    # this function is deprecated in more_itertools
    # because it is superseded by the standard library
    from more_itertools import zip_equal

    def zip(*args, strict=False):
        if strict:
            return zip_equal(*args)
        else:
            return builtin_zip(*args)

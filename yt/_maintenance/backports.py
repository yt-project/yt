# A namespace to store simple backports from most recent Python versions
# Backports should be contained in (if/else) blocks, cheking on the runtime
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

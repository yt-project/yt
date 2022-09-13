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

if sys.version_info < (3, 8):
    from _thread import RLock  # type: ignore
    from typing import List

    GenericAlias = type(List[int])

    _NOT_FOUND = object()

    class cached_property:
        def __init__(self, func):
            self.func = func
            self.attrname = None
            self.__doc__ = func.__doc__
            self.lock = RLock()

        def __set_name__(self, owner, name):
            if self.attrname is None:
                self.attrname = name
            elif name != self.attrname:
                raise TypeError(
                    "Cannot assign the same cached_property to two different names "
                    f"({self.attrname!r} and {name!r})."
                )

        def __get__(self, instance, owner=None):
            if instance is None:
                return self
            if self.attrname is None:
                raise TypeError(
                    "Cannot use cached_property instance without calling __set_name__ on it."
                )
            try:
                cache = instance.__dict__
            except AttributeError:  # not all objects have __dict__ (e.g. class defines slots)
                msg = (
                    f"No '__dict__' attribute on {type(instance).__name__!r} "
                    f"instance to cache {self.attrname!r} property."
                )
                raise TypeError(msg) from None
            val = cache.get(self.attrname, _NOT_FOUND)
            if val is _NOT_FOUND:
                with self.lock:
                    # check if another thread filled cache while we awaited lock
                    val = cache.get(self.attrname, _NOT_FOUND)
                    if val is _NOT_FOUND:
                        val = self.func(instance)
                        try:
                            cache[self.attrname] = val
                        except TypeError:
                            msg = (
                                f"The '__dict__' attribute on {type(instance).__name__!r} instance "
                                f"does not support item assignment for caching {self.attrname!r} property."
                            )
                            raise TypeError(msg) from None
            return val

        __class_getitem__ = classmethod(GenericAlias)

else:
    pass


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

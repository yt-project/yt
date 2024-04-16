# avoid deprecation warnings in numpy >= 2.0

from importlib.metadata import version

from packaging.version import Version

NUMPY_VERSION = Version(version("numpy"))

if NUMPY_VERSION >= Version("2.0.0dev0"):
    from numpy import trapezoid as trapezoid  # type: ignore [attr-defined]
else:
    from numpy import trapz as trapezoid  # noqa: F401

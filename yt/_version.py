from typing import NamedTuple

from packaging.version import Version

__all__ = [
    "__version__",
    "version_info",
]

__version__ = "4.1.4"  # keep in sync with setup.cfg


class VersionTuple(NamedTuple):
    """
    A minimal representation of the current version number
    that can be used downstream to check the runtime version
    simply by comparing with builtin tuples, as can be done with
    the runtime Python version using sys.version_info

    https://docs.python.org/3/library/sys.html#sys.version_info
    """

    major: int
    minor: int
    micro: int
    releaselevel: str
    serial: int


def _parse_to_version_info(version_str: str) -> VersionTuple:
    # adapted from matplotlib 3.5
    """
    Parse a version string to a namedtuple analogous to sys.version_info.
    See:
    https://packaging.pypa.io/en/latest/version.html#packaging.version.parse
    https://docs.python.org/3/library/sys.html#sys.version_info
    """
    v = Version(version_str)
    if v.pre is None and v.post is None and v.dev is None:
        return VersionTuple(v.major, v.minor, v.micro, "final", 0)
    elif v.dev is not None:
        return VersionTuple(v.major, v.minor, v.micro, "alpha", v.dev)
    elif v.pre is not None:
        releaselevel = {"a": "alpha", "b": "beta", "rc": "candidate"}.get(
            v.pre[0], "alpha"
        )
        return VersionTuple(v.major, v.minor, v.micro, releaselevel, v.pre[1])
    elif v.post is not None:
        # fallback for v.post: guess-next-dev scheme from setuptools_scm
        return VersionTuple(v.major, v.minor, v.micro + 1, "alpha", v.post)
    else:
        return VersionTuple(v.major, v.minor, v.micro + 1, "alpha", 0)


version_info = _parse_to_version_info(__version__)

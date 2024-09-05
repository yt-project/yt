import pytest

import yt
from yt._version import VersionTuple, _parse_to_version_info


@pytest.mark.parametrize(
    "version_str, expected",
    (
        ("4.1.0", VersionTuple(4, 1, 0, "final", 0)),
        ("6.2.5", VersionTuple(6, 2, 5, "final", 0)),
        ("4.1.dev0", VersionTuple(4, 1, 0, "alpha", 0)),
        ("4.1.0rc", VersionTuple(4, 1, 0, "candidate", 0)),
        ("4.1.0rc1", VersionTuple(4, 1, 0, "candidate", 1)),
        ("4.1.0rc2", VersionTuple(4, 1, 0, "candidate", 2)),
    ),
)
def test_parse_version_info(version_str, expected):
    actual = _parse_to_version_info(version_str)
    assert actual == expected


def test_version_tuple_comp():
    # exercise comparison with builtin tuples
    # comparison results do not matter for this test

    yt.version_info > (4,)  # noqa: B015
    yt.version_info > (4, 1)  # noqa: B015
    yt.version_info > (4, 1, 0)  # noqa: B015
    yt.version_info < (4, 1, 0)  # noqa: B015
    yt.version_info <= (4, 1, 0)  # noqa: B015
    yt.version_info >= (4, 1, 0)  # noqa: B015
    yt.version_info == (4, 1, 0)  # noqa: B015

    assert isinstance(yt.version_info, tuple)

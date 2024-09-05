import re

import pytest

from yt.testing import fake_amr_ds
from yt.visualization.plot_window import SlicePlot


@pytest.mark.parametrize(
    ("origin", "msg"),
    [
        (
            "ONE",
            (
                "Invalid origin argument. "
                "Single element specification must be 'window', 'domain', or 'native'. "
                "Got 'ONE'"
            ),
        ),
        (
            "ONE-TWO",
            (
                "Invalid origin argument. "
                "Using 2 elements:\n"
                " - the first one must be 'left', 'right', 'lower', 'upper' or 'center'\n"
                " - the second one must be 'window', 'domain', or 'native'\n"
                "Got 'ONE-TWO'"
            ),
        ),
        (
            "ONE-window",
            (
                "Invalid origin argument. "
                "Using 2 elements:\n"
                " - the first one must be 'left', 'right', 'lower', 'upper' or 'center'\n"
                "Got 'ONE-window'"
            ),
        ),
        (
            "left-TWO",
            (
                "Invalid origin argument. "
                "Using 2 elements:\n"
                " - the second one must be 'window', 'domain', or 'native'\n"
                "Got 'left-TWO'"
            ),
        ),
        (
            "ONE-TWO-THREE",
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                " - the second one must be 'left', 'right', 'center' or a distance\n"
                " - the third one must be 'window', 'domain', or 'native'\n"
                "Got 'ONE-TWO-THREE'"
            ),
        ),
        (
            "ONE-TWO-window",
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                " - the second one must be 'left', 'right', 'center' or a distance\n"
                "Got 'ONE-TWO-window'"
            ),
        ),
        (
            "ONE-left-window",
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                "Got 'ONE-left-window'"
            ),
        ),
        (
            "ONE-left-THREE",
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                " - the third one must be 'window', 'domain', or 'native'\n"
                "Got 'ONE-left-THREE'"
            ),
        ),
        (
            "lower-left-THREE",
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the third one must be 'window', 'domain', or 'native'\n"
                "Got 'lower-left-THREE'"
            ),
        ),
        (
            ("ONE", "TWO", "THREE"),
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                " - the second one must be 'left', 'right', 'center' or a distance\n"
                " - the third one must be 'window', 'domain', or 'native'\n"
                "Got ('ONE', 'TWO', 'THREE')"
            ),
        ),
        (
            ("ONE", "TWO", (1, 1, 3)),
            (
                "Invalid origin argument. "
                "Using 3 elements:\n"
                " - the first one must be 'lower', 'upper' or 'center' or a distance\n"
                " - the second one must be 'left', 'right', 'center' or a distance\n"
                " - the third one must be 'window', 'domain', or 'native'\n"
                "Got ('ONE', 'TWO', (1, 1, 3))"
            ),
        ),
        (
            "ONE-TWO-THREE-FOUR",
            (
                "Invalid origin argument with too many elements; "
                "expected 1, 2 or 3 elements, got 'ONE-TWO-THREE-FOUR', counting 4 elements. "
                "Use '-' as a separator for string arguments."
            ),
        ),
    ],
)
def test_invalidate_origin_value(origin, msg):
    ds = fake_amr_ds(fields=[("gas", "density")], units=["g*cm**-3"])
    with pytest.raises(ValueError, match=re.escape(msg)):
        SlicePlot(ds, "z", ("gas", "density"), origin=origin)


@pytest.mark.parametrize(
    "origin",
    # don't attempt to match exactly a TypeError message because it should be
    # emitted by unyt, not yt
    [
        ((50, 50, 50), "TWO", "THREE"),
        ("ONE", (50, 50, 50), "THREE"),
    ],
)
def test_invalidate_origin_type(origin):
    ds = fake_amr_ds(fields=[("gas", "density")], units=["g*cm**-3"])
    with pytest.raises(TypeError):
        SlicePlot(ds, "z", ("gas", "density"), origin=origin)

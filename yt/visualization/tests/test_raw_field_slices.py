import yt
from yt.utilities.answer_testing.framework import (
    GenericImageTest,
    data_dir_load,
    requires_ds,
)


def setup():
    """Test specific setup."""
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def compare(ds, field, test_prefix, decimals=12):
    def slice_image(filename_prefix):
        sl = yt.SlicePlot(ds, "z", field)
        sl.set_log("all", False)
        image_file = sl.save(filename_prefix)
        return image_file

    slice_image.__name__ = f"slice_{test_prefix}"
    test = GenericImageTest(ds, slice_image, decimals)
    test.prefix = test_prefix
    return test


raw_fields = "Laser/plt00015"
_raw_field_names = [
    ("raw", "Bx"),
    ("raw", "By"),
    ("raw", "Bz"),
    ("raw", "Ex"),
    ("raw", "Ey"),
    ("raw", "Ez"),
    ("raw", "jx"),
    ("raw", "jy"),
    ("raw", "jz"),
]


@requires_ds(raw_fields)
def test_raw_field_slices():
    ds = data_dir_load(raw_fields)
    for field in _raw_field_names:
        yield compare(ds, field, f"answers_raw_{field[1]}")

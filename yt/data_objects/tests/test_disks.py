import pytest

from yt import YTQuantity
from yt.testing import fake_random_ds


def test_bad_disk_input():
    # Fixes 1768
    ds = fake_random_ds(16)

    # Test invalid 3d array
    with pytest.raises(
        TypeError,
        match=r"^Expected an array of size \(3,\), received 'list' of length 4$",
    ):
        ds.disk(ds.domain_center, [0, 0, 1, 1], (10, "kpc"), (20, "kpc"))

    # Test invalid float
    with pytest.raises(
        TypeError,
        match=(
            r"^Expected a numeric value \(or size-1 array\), "
            r"received 'unyt.array.unyt_array' of length 3$"
        ),
    ):
        ds.disk(ds.domain_center, [0, 0, 1], ds.domain_center, (20, "kpc"))

    # Test invalid float
    with pytest.raises(
        TypeError,
        match=(
            r"^Expected a numeric value \(or tuple of format \(float, String\)\), "
            r"received an inconsistent tuple '\(10, 10\)'.$"
        ),
    ):
        ds.disk(ds.domain_center, [0, 0, 1], (10, 10), (20, "kpc"))

    # Test invalid iterable
    with pytest.raises(
        TypeError,
        match=r"^Expected an iterable object, received 'unyt\.array\.unyt_quantity'$",
    ):
        ds.disk(
            ds.domain_center,
            [0, 0, 1],
            (10, "kpc"),
            (20, "kpc"),
            fields=YTQuantity(1, "kpc"),
        )

    # Test invalid object
    with pytest.raises(
        TypeError,
        match=(
            r"^Expected an object of 'yt\.data_objects\.static_output\.Dataset' type, "
            r"received 'yt\.data_objects\.selection_objects\.region\.YTRegion'$"
        ),
    ):
        ds.disk(ds.domain_center, [0, 0, 1], (10, "kpc"), (20, "kpc"), ds=ds.all_data())

    # Test valid disk
    ds.disk(ds.domain_center, [0, 0, 1], (10, "kpc"), (20, "kpc"))
    ds.disk(ds.domain_center, [0, 0, 1], 10, (20, "kpc"))

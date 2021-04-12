from collections import defaultdict
from itertools import product

import pytest

from yt.testing import fake_random_ds
from yt.utilities.exceptions import YTFieldNotFound


def test_unexisting_field_access():
    ds = fake_random_ds(16, particles=10)

    fname2ftype = defaultdict(list)

    for ft, fn in ds.derived_field_list:
        fname2ftype[fn].append(ft)

    ad = ds.all_data()

    ftypes = ("gas", "io")
    fnames = (
        "density",
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
    )

    # Try invalid ftypes, fnames combinations
    for ft, fn in product(ftypes, fnames):
        if (ft, fn) in ds.derived_field_list:
            continue

        with pytest.raises(YTFieldNotFound) as excinfo:
            ad[(ft, fn)]

        # Make sure the existing field has been suggested
        for possible_ft in fname2ftype[fn]:
            assert str((possible_ft, fn)) in str(excinfo.value)

    # Try typos
    for bad_field, good_field in (
        (("gas", "densi_y"), ("gas", "density")),
        (("oi", "particle_mass"), ("io", "particle_mass")),
        (("gas", "DENSITY"), ("gas", "density")),
    ):
        with pytest.raises(YTFieldNotFound) as excinfo:
            ad[bad_field]

        assert str(good_field) in str(excinfo.value)

from collections import defaultdict

import pytest

from yt._maintenance.deprecation import VisibleDeprecationWarning
from yt.testing import fake_amr_ds


def test_ambiguous_fails():
    ds = fake_amr_ds(particles=10)
    msg = "The requested field name '{}' is ambiguous"

    fnames_to_ftype = defaultdict(list)
    for ftype, fname in ds.field_list:
        fnames_to_ftype[fname].append(ftype)

    ambiguous_fnames = [
        fname for fname, ftypes in fnames_to_ftype.items() if len(ftypes) > 1
    ]
    unambiguous_fnames = [
        fname for fname, ftypes in fnames_to_ftype.items() if len(ftypes) == 1
    ] + ["zeros", "ones"]

    # Make sure we are testing both cases with a few cases
    assert len(ambiguous_fnames) > 2
    assert len(unambiguous_fnames) > 2

    # Test warnings are issued for ambiguous fields
    for fname in ambiguous_fnames:
        with pytest.warns(VisibleDeprecationWarning, match=msg.format(fname)):
            ds.r[fname]

    # Test no warnings are issued for single fname access that aren't ambiguous
    for fname in unambiguous_fnames:
        with pytest.warns(None) as record:
            ds.r[fname]
        assert len(record) == 0

    # Test no warning are issued for tuple access
    for ftype, fname in ds.field_list:
        with pytest.warns(None) as record:
            ds.r[(ftype, fname)]

        assert len(record) == 0

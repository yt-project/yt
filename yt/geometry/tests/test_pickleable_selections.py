import pickle

import numpy as np

from yt.testing import assert_equal, fake_particle_ds


def test_pickleability():
    # tests the pickleability of the selection objects.

    def pickle_test(sel_obj):
        assert_fields = sel_obj._get_state_attnames()  # the attrs used in get/set state
        new_sel_obj = pickle.loads(pickle.dumps(sel_obj))
        for attr in assert_fields:
            assert_equal(getattr(new_sel_obj, attr), getattr(sel_obj, attr))

    # list of selection types and argument tuples for each selection type
    c = np.array([0.5, 0.5, 0.5])
    sargs = (
        ("point", (c,)),
        ("sphere", (c, 0.25)),
        ("box", (c - 0.3, c)),
        ("ellipsoid", (c, 0.3, 0.2, 0.1, c - 0.4, 0.2)),
        ("disk", (c, [1, 0, 0], 0.2, 0.2)),
        ("cutting", ([0.1, 0.2, -0.9], [0.5, 0.42, 0.6])),
        ("ortho_ray", ("z", c)),
        ("ray", (c, [0.1, 0.1, 0.1])),
    )

    # load fake data
    ds = fake_particle_ds()
    for sel_type, args in sargs:
        sel = getattr(ds, sel_type)(*args)  # instantiate this selection type
        pickle_test(sel.selector)  # make sure it (un)pickles

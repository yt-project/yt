import pickle

import numpy as np

from yt.testing import (  # fake_amr_ds,; fake_random_ds,; requires_module,
    assert_equal,
    fake_particle_ds,
)


def test_pickleability():
    # only tests the pickleability of the selection objects.

    def pickle_test(sel_obj, assert_fields=[]):
        new_sel_obj = pickle.loads(pickle.dumps(sel_obj))
        for attr in assert_fields:
            assert_equal(getattr(new_sel_obj, attr), getattr(sel_obj, attr))

    ds = fake_particle_ds()
    c = [0.5, 0.5, 0.5]
    pickle_test(ds.point(c).selector, ["p"])
    pickle_test(ds.sphere(c, 0.25).selector, ["center", "radius", "radius2"])
    sel_ob = ds.box([0.2, 0.2, 0.2], c)
    pickle_test(sel_ob.selector, ["left_edge", "right_edge", "is_all_data"])
    sel_ob = ds.ellipsoid(c, 0.3, 0.2, 0.1, np.array([0.1, 0.1, 0.1]), 0.2)
    pickle_test(sel_ob.selector, ["vec", "center", "mag"])
    sel_ob = ds.disk(c, [1, 0, 0], 0.2, 0.2)
    pickle_test(sel_ob.selector, ["center", "radius", "radius2", "norm_vec", "height"])
    sel_ob = ds.cutting([0.1, 0.2, -0.9], [0.5, 0.42, 0.6])
    pickle_test(sel_ob.selector, ["d", "norm_vec"])


# def test_unpickled_selections():
#     # TO DO: use fake dataset to test selections
#     return

# TO DO: add some answer test

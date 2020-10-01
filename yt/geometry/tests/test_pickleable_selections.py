import pickle

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
    pickle_test(
        ds.box([0.2, 0.2, 0.2], c).selector, ["left_edge", "right_edge", "is_all_data"]
    )


# def test_unpickled_selections():
#     # TO DO: use fake dataset to test selections
#     return

# TO DO: add some answer test

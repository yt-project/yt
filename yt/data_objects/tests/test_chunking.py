from yt.testing import assert_equal, assert_true, fake_random_ds
from yt.units._numpy_wrapper_functions import uconcatenate


def _get_dobjs(c):
    dobjs = [
        ("sphere", ("center", (1.0, "unitary"))),
        ("sphere", ("center", (0.1, "unitary"))),
        ("ortho_ray", (0, (c[1], c[2]))),
        ("slice", (0, c[0])),
        # ("disk", ("center", [0.1, 0.3, 0.6],
        #           (0.2, 'unitary'), (0.1, 'unitary'))),
        ("cutting", ([0.1, 0.3, 0.6], "center")),
        ("all_data", ()),
    ]
    return dobjs


def test_chunking():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs=nprocs)
        c = (ds.domain_right_edge + ds.domain_left_edge) / 2.0
        c += ds.arr(0.5 / ds.domain_dimensions, "code_length")
        for dobj in _get_dobjs(c):
            obj = getattr(ds, dobj[0])(*dobj[1])
            coords = {"f": {}, "i": {}}
            for t in ["io", "all", "spatial"]:
                coords["i"][t] = []
                coords["f"][t] = []
                for chunk in obj.chunks(None, t):
                    coords["f"][t].append(chunk.fcoords[:, :])
                    coords["i"][t].append(chunk.icoords[:, :])
                coords["f"][t] = uconcatenate(coords["f"][t])
                coords["i"][t] = uconcatenate(coords["i"][t])
                coords["f"][t].sort()
                coords["i"][t].sort()
            assert_equal(coords["f"]["io"], coords["f"]["all"])
            assert_equal(coords["f"]["io"], coords["f"]["spatial"])
            assert_equal(coords["i"]["io"], coords["i"]["all"])
            assert_equal(coords["i"]["io"], coords["i"]["spatial"])


def test_ds_hold():
    ds1 = fake_random_ds(64)
    ds2 = fake_random_ds(128)
    dd = ds1.all_data()
    assert_true(
        dd.ds.__hash__() == ds1.__hash__()
    )  # dd.ds is a weakref, so can't use "is"
    assert_true(dd.index is ds1.index)
    assert_equal(dd[("index", "ones")].size, 64**3)
    with dd._ds_hold(ds2):
        assert_true(dd.ds.__hash__() == ds2.__hash__())
        assert_true(dd.index is ds2.index)
        assert_equal(dd[("index", "ones")].size, 128**3)
    assert_true(dd.ds.__hash__() == ds1.__hash__())
    assert_true(dd.index is ds1.index)
    assert_equal(dd[("index", "ones")].size, 64**3)

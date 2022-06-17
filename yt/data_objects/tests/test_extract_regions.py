import numpy as np

from yt.loaders import load
from yt.testing import (
    assert_almost_equal,
    assert_equal,
    fake_amr_ds,
    fake_random_ds,
    requires_file,
)


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


def test_cut_region():
    # We decompose in different ways
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(
            64,
            nprocs=nprocs,
            fields=("density", "temperature", "velocity_x"),
            units=("g/cm**3", "K", "cm/s"),
        )
        # We'll test two objects
        dd = ds.all_data()
        r = dd.cut_region(
            [
                "obj[('gas', 'temperature')] > 0.5",
                "obj[('gas', 'density')] < 0.75",
                "obj[('gas', 'velocity_x')] > 0.25",
            ]
        )
        t = (
            (dd[("gas", "temperature")] > 0.5)
            & (dd[("gas", "density")] < 0.75)
            & (dd[("gas", "velocity_x")] > 0.25)
        )
        assert_equal(np.all(r[("gas", "temperature")] > 0.5), True)
        assert_equal(np.all(r[("gas", "density")] < 0.75), True)
        assert_equal(np.all(r[("gas", "velocity_x")] > 0.25), True)
        assert_equal(np.sort(dd[("gas", "density")][t]), np.sort(r[("gas", "density")]))
        assert_equal(np.sort(dd[("index", "x")][t]), np.sort(r[("index", "x")]))
        r2 = r.cut_region(["obj[('gas', 'temperature')] < 0.75"])
        t2 = r[("gas", "temperature")] < 0.75
        assert_equal(
            np.sort(r2[("gas", "temperature")]), np.sort(r[("gas", "temperature")][t2])
        )
        assert_equal(np.all(r2[("gas", "temperature")] < 0.75), True)

        # Now we can test some projections
        dd = ds.all_data()
        cr = dd.cut_region(["obj[('index', 'ones')] > 0"])
        for weight in [None, ("gas", "density")]:
            p1 = ds.proj(("gas", "density"), 0, data_source=dd, weight_field=weight)
            p2 = ds.proj(("gas", "density"), 0, data_source=cr, weight_field=weight)
            for f in p1.field_data:
                assert_almost_equal(p1[f], p2[f])
        cr = dd.cut_region(["obj[('gas', 'density')] > 0.25"])
        p2 = ds.proj(("gas", "density"), 2, data_source=cr)
        assert_equal(p2[("gas", "density")].max() > 0.25, True)
        p2 = ds.proj(
            ("gas", "density"), 2, data_source=cr, weight_field=("gas", "density")
        )
        assert_equal(p2[("gas", "density")].max() > 0.25, True)


def test_region_and_particles():
    ds = fake_amr_ds(particles=10000)

    ad = ds.all_data()
    reg = ad.cut_region('obj[("index", "x")] < .5')

    mask = ad[("all", "particle_position_x")] < 0.5
    expected = np.sort(ad[("all", "particle_position_x")][mask].value)
    result = np.sort(reg[("all", "particle_position_x")])

    assert_equal(expected.shape, result.shape)
    assert_equal(expected, result)


ISOGAL = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(ISOGAL)
def test_region_chunked_read():
    # see #2104
    ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")

    sp = ds.sphere((0.5, 0.5, 0.5), (2, "kpc"))
    dense_sp = sp.cut_region(['obj[("gas", "H_p0_number_density")]>= 1e-2'])
    dense_sp.quantities.angular_momentum_vector()


@requires_file(ISOGAL)
def test_chained_cut_region():
    # see Issue #2233
    ds = load(ISOGAL)
    base = ds.disk([0.5, 0.5, 0.5], [0, 0, 1], (4, "kpc"), (10, "kpc"))
    c1 = "(obj[('index', 'cylindrical_radius')].in_units('kpc') > 2.0)"
    c2 = "(obj[('gas', 'density')].to('g/cm**3') > 1e-26)"

    cr12 = base.cut_region([c1, c2])
    cr1 = base.cut_region([c1])
    cr12c = cr1.cut_region([c2])

    field = ("index", "cell_volume")
    assert_equal(
        cr12.quantities.total_quantity(field), cr12c.quantities.total_quantity(field)
    )

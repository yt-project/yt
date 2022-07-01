from yt.testing import assert_allclose_units, fake_amr_ds


def test_icoords_to_ires():
    for geometry in ("cartesian", "spherical", "cylindrical"):
        ds = fake_amr_ds(geometry=geometry)

        dd = ds.r[:]
        icoords = dd.icoords
        ires = dd.ires
        fcoords, fwidth = ds.index._icoords_to_fcoords(icoords, ires)
        assert_allclose_units(fcoords, dd.fcoords, rtol=1e-14)
        assert_allclose_units(fwidth, dd.fwidth, rtol=1e-14)

        fcoords_xz, fwidth_xz = ds.index._icoords_to_fcoords(
            icoords[:, (0, 2)], ires, axes=(0, 2)
        )
        assert_allclose_units(fcoords_xz[:, 0], dd.fcoords[:, 0])
        assert_allclose_units(fcoords_xz[:, 1], dd.fcoords[:, 2])
        assert_allclose_units(fwidth_xz[:, 0], dd.fwidth[:, 0])
        assert_allclose_units(fwidth_xz[:, 1], dd.fwidth[:, 2])

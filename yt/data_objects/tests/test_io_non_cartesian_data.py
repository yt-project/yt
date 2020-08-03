from tempfile import NamedTemporaryFile

from yt.frontends.ytdata.data_structures import YTDataContainerDataset
from yt.testing import fake_amr_ds


def test_preserve_geometric_properties():
    for geom in ("cartesian", "cylindrical", "spherical"):
        ds1 = fake_amr_ds(fields=[("gas", "density")], geometry=geom)
        ad = ds1.all_data()
        with NamedTemporaryFile(suffix=".h5") as tmpf:
            fn = ad.save_as_dataset(tmpf.name, fields=["density"])
            ds2 = YTDataContainerDataset(fn)
            dfl = ds2.derived_field_list
        assert ds1.geometry == ds2.geometry == geom

        expected = set(ds1.coordinates.axis_order)
        actual = set([fname for ftype, fname in dfl])
        assert expected.difference(actual) == set()

import os
from tempfile import TemporaryDirectory

import numpy as np

from yt.frontends.ytdata.api import save_as_dataset
from yt.frontends.ytdata.data_structures import YTDataContainerDataset
from yt.loaders import load
from yt.testing import fake_amr_ds, requires_module
from yt.units import YTQuantity


@requires_module("h5py")
def test_preserve_geometric_properties():
    for geom in ("cartesian", "cylindrical", "spherical"):
        ds1 = fake_amr_ds(fields=[("gas", "density")], units=["g/cm**3"], geometry=geom)
        ad = ds1.all_data()
        with TemporaryDirectory() as tmpdir:
            tmpf = os.path.join(tmpdir, "savefile.h5")
            fn = ad.save_as_dataset(tmpf, fields=[("gas", "density")])
            ds2 = load(fn)
            assert isinstance(ds2, YTDataContainerDataset)
            dfl = ds2.derived_field_list
        assert ds1.geometry == ds2.geometry == geom

        expected = set(ds1.coordinates.axis_order)
        actual = {fname for ftype, fname in dfl}
        assert expected.difference(actual) == set()


@requires_module("h5py")
def test_default_to_cartesian():
    data = {"density": np.random.random(128)}
    ds_attrs = {"current_time": YTQuantity(10, "Myr")}
    with TemporaryDirectory() as tmpdir:
        tmpf = os.path.join(tmpdir, "savefile.h5")
        fn = save_as_dataset(ds_attrs, tmpf, data)
        ds2 = load(fn)
    assert ds2.geometry == "cartesian"

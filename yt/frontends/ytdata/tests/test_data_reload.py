import numpy as np

from yt.loaders import load_uniform_grid
from yt.testing import requires_module_pytest
from yt.utilities.on_demand_imports import _h5py as h5py


@requires_module_pytest("h5py")
def test_save_as_data_unit_system(tmp_path):
    # This test checks that the file saved with calling save_as_dataset
    # for a ds with a "code" unit system contains the proper "unit_system_name".
    # It checks the hdf5 file directly rather than using yt.load(), because
    # https://github.com/yt-project/yt/issues/4315 only manifested restarting
    # the python kernel (because the unit registry is state dependent).

    fi = tmp_path / "output_data.h5"
    shp = (4, 4, 4)
    data = {"density": np.random.random(shp)}
    ds = load_uniform_grid(data, shp, unit_system="code")
    assert "code" in ds._unit_system_name

    sp = ds.sphere(ds.domain_center, ds.domain_width[0] / 2.0)
    sp.save_as_dataset(fi)

    with h5py.File(fi, mode="r") as f:
        assert f.attrs["unit_system_name"] == "code"

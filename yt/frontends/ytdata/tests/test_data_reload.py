import numpy as np

from yt.loaders import load, load_uniform_grid
from yt.testing import fake_amr_ds, requires_module_pytest
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


@requires_module_pytest("h5py")
def test_save_as_data_chunk(tmp_path):
    """
    Tests whether a saved (sphere) dataset matches what was saved, and whether
    accessing the saved data via different paths produces the same result.
    """
    sphere_path = tmp_path / "test_sphere.h5"
    ds = fake_amr_ds(length_unit="kpc")
    sp = ds.sphere(ds.domain_center, (0.25, "kpc"))
    original_data = sp["stream", "Density"]
    sp.save_as_dataset(sphere_path, fields=[("stream", "Density")])

    sp_ds = load(sphere_path)  # should produce a 7-chunk dataset
    assert len(sp_ds.index.data_files) > 1, "Test data not chunked."

    reloaded_data = sp_ds.data["grid", "Density"]  # previously crashed with chunking
    all_reloaded_data = sp_ds.all_data()["grid", "Density"]

    np.testing.assert_array_equal(
        original_data,
        reloaded_data,
        err_msg="Reloaded sphere produces different data to what was saved.",
    )

    np.testing.assert_array_equal(
        all_reloaded_data,
        reloaded_data,
        err_msg="Reloaded sphere produces different results from .data vs .all_data().",
    )

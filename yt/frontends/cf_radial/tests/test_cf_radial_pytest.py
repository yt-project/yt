from pathlib import Path

import numpy as np

import yt
from yt.frontends.cf_radial.api import CFRadialDataset
from yt.testing import requires_module_pytest as requires_module


def create_cf_radial_mock_gridded_ds(savedir: Path) -> Path:
    # a minimal dataset that yt should pick up as a CfRadial dataset
    import xarray as xr

    file_to_save = savedir / "mock_gridded_cfradial.nc"

    shp = (16, 16, 16)
    xyz = {dim: np.linspace(0.0, 1.0, shp[idim]) for idim, dim in enumerate("xyz")}
    da = xr.DataArray(data=np.ones(shp), coords=xyz, name="reflectivity")
    ds_xr = da.to_dataset()
    ds_xr.attrs["conventions"] = "CF/Radial"
    ds_xr.x.attrs["units"] = "m"
    ds_xr.y.attrs["units"] = "m"
    ds_xr.z.attrs["units"] = "m"
    ds_xr.reflectivity.attrs["units"] = ""

    times = np.array(["2017-05-19T01", "2017-05-19T01:01"], dtype="datetime64[ns]")
    ds_xr = ds_xr.assign_coords({"time": times})
    ds_xr["origin_latitude"] = xr.DataArray(np.zeros_like(times), dims=("time",))
    ds_xr["origin_longitude"] = xr.DataArray(np.zeros_like(times), dims=("time",))

    ds_xr.to_netcdf(file_to_save, engine="netcdf4")

    return file_to_save


@requires_module("xarray", "netCDF4", "pyart")
def test_load_mock_gridded_cf_radial(tmp_path):
    import xarray as xr

    test_file = create_cf_radial_mock_gridded_ds(tmp_path)
    assert test_file.exists()

    # make sure that the mock dataset is valid and can be re-loaded with xarray
    with xr.open_dataset(test_file) as ds_xr:
        assert "CF/Radial" in ds_xr.conventions

    test_file = create_cf_radial_mock_gridded_ds(tmp_path)

    ds_yt = yt.load(test_file)
    assert isinstance(ds_yt, CFRadialDataset)

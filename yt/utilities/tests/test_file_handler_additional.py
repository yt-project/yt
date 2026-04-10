import os
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from numpy.testing import assert_raises

import yt.utilities.file_handler as fh
import yt.utilities.on_demand_imports as odi


def _write_bytes(path, payload):
    with open(path, "wb") as stream:
        stream.write(payload)


class _FakeHandle(dict):
    def __init__(self):
        super().__init__({"group": 3.0})
        self.attrs = {"units": "K"}
        self.closed = False

    def close(self):
        self.closed = True


class _FakeImageBaseHDU:
    pass


class _FakeHDUList(list):
    def __init__(self, value=None):
        if value is None:
            super().__init__()
        elif isinstance(value, list):
            super().__init__(value)
        else:
            super().__init__([value])
        self.closed = False

    def close(self):
        self.closed = True


class _FakeClosable:
    def __init__(self):
        self.closed = False

    def close(self):
        self.closed = True


class _FakeNetCDFDataset:
    def __init__(self, filename, mode, kwargs):
        self.filename = filename
        self.mode = mode
        self.kwargs = kwargs
        self.closed = False

    def close(self):
        self.closed = True


def test_file_signatures_and_deprecated_warning_guards():
    with TemporaryDirectory() as tmpdir:
        hdf5_path = os.path.join(tmpdir, "sample.h5")
        classic_nc_path = os.path.join(tmpdir, "sample.nc")
        netcdf4_path = os.path.join(tmpdir, "sample.nc4")
        text_path = os.path.join(tmpdir, "sample.txt")
        missing_path = os.path.join(tmpdir, "missing.dat")

        _write_bytes(hdf5_path, b"\x89HDF\r\n\x1a\npayload")
        _write_bytes(classic_nc_path, b"CDF\x01payload")
        _write_bytes(netcdf4_path, b"\x89HDF\r\n\x1a\npayload")
        _write_bytes(text_path, b"TEXTpayload")

        assert fh.valid_hdf5_signature(hdf5_path) is True
        assert fh.valid_hdf5_signature(text_path) is False
        assert fh.valid_hdf5_signature(missing_path) is False

        assert fh.valid_netcdf_signature(classic_nc_path) is True
        assert fh.valid_netcdf_signature(netcdf4_path) is True
        assert fh.valid_netcdf_signature(text_path) is False
        assert fh.valid_netcdf_signature(missing_path) is False

        original_warning = fh.issue_deprecation_warning
        original_h5py = fh.h5py
        original_netcdf4 = odi._netCDF4
        try:
            fh.issue_deprecation_warning = lambda *args, **kwargs: None
            fh.h5py = SimpleNamespace(File=fh.NotAModule("h5py"))
            odi._netCDF4 = SimpleNamespace(Dataset=fh.NotAModule("netCDF4"))

            assert fh.valid_netcdf_classic_signature(classic_nc_path) is True
            assert fh.valid_netcdf_classic_signature(text_path) is False
            assert fh.valid_netcdf_classic_signature(missing_path) is False

            with assert_raises(RuntimeError):
                fh.warn_h5py(hdf5_path)
            fh.warn_h5py(text_path)

            with assert_raises(RuntimeError):
                fh.warn_netcdf(netcdf4_path)
            fh.warn_netcdf(text_path)
        finally:
            fh.issue_deprecation_warning = original_warning
            fh.h5py = original_h5py
            odi._netCDF4 = original_netcdf4


def test_hdf5_file_handler_wraps_underlying_handle_methods():
    fake_handle = _FakeHandle()
    original_h5py = fh.h5py
    try:
        fh.h5py = SimpleNamespace(File=lambda filename, mode="r": fake_handle)

        handler = fh.HDF5FileHandler("dummy.h5")

        assert handler["group"] == 3.0
        assert "group" in handler
        assert len(handler) == 1
        assert handler.attrs == {"units": "K"}
        assert handler.keys() == ["group"]
        assert handler.items() == [("group", 3.0)]

        handler.close()
        assert fake_handle.closed is True

        handler.handle = None
        handler.close()
    finally:
        fh.h5py = original_h5py


def test_fits_file_handler_covers_all_constructor_branches_and_cleanup():
    open_calls = []

    def fake_open(filename, **kwargs):
        open_calls.append((filename, kwargs))
        return _FakeHDUList(["opened"])

    fake_astropy = SimpleNamespace(
        pyfits=SimpleNamespace(
            hdu=SimpleNamespace(image=SimpleNamespace(_ImageBaseHDU=_FakeImageBaseHDU)),
            HDUList=_FakeHDUList,
            open=fake_open,
        )
    )

    original_astropy = odi._astropy
    try:
        odi._astropy = fake_astropy

        image_hdu = _FakeImageBaseHDU()
        image_handler = fh.FITSFileHandler(image_hdu)
        assert isinstance(image_handler.handle, _FakeHDUList)
        image_handler.close()
        assert image_handler.handle.closed is True

        hdu_list = _FakeHDUList(["existing"])
        list_handler = fh.FITSFileHandler(hdu_list)
        assert list_handler.handle is hdu_list

        path_handler = fh.FITSFileHandler("image.fits")
        assert open_calls == [
            (
                "image.fits",
                {
                    "memmap": True,
                    "do_not_scale_image_data": True,
                    "ignore_blank": True,
                },
            )
        ]

        extra_a = _FakeClosable()
        extra_b = _FakeClosable()
        path_handler._fits_files = [extra_a, extra_b]
        path_handler.__del__()
        path_handler._fits_files = []
        assert extra_a.closed is True
        assert extra_b.closed is True
        assert path_handler.handle is None
    finally:
        odi._astropy = original_astropy


def test_netcdf4_file_handler_context_manager_closes_dataset():
    created = []

    def fake_dataset(filename, mode="r", **kwargs):
        dataset = _FakeNetCDFDataset(filename, mode, kwargs)
        created.append(dataset)
        return dataset

    original_netcdf4 = odi._netCDF4
    try:
        odi._netCDF4 = SimpleNamespace(Dataset=fake_dataset)

        handler = fh.NetCDF4FileHandler("sample.nc")
        with handler.open_ds(maskandscale=True) as dataset:
            assert dataset.filename == "sample.nc"
            assert dataset.mode == "r"
            assert dataset.kwargs == {"maskandscale": True}

        assert len(created) == 1
        assert created[0].closed is True
    finally:
        odi._netCDF4 = original_netcdf4

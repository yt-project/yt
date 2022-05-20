from contextlib import contextmanager

from yt.utilities.on_demand_imports import NotAModule, _h5py as h5py


def valid_hdf5_signature(fn):
    signature = b"\x89HDF\r\n\x1a\n"
    try:
        with open(fn, "rb") as f:
            header = f.read(8)
            return header == signature
    except Exception:
        return False


def warn_h5py(fn):
    needs_h5py = valid_hdf5_signature(fn)
    if needs_h5py and isinstance(h5py.File, NotAModule):
        raise RuntimeError(
            "This appears to be an HDF5 file, but h5py is not installed."
        )


class HDF5FileHandler:
    handle = None

    def __init__(self, filename):
        self.handle = h5py.File(filename, mode="r")

    def __getitem__(self, key):
        return self.handle[key]

    def __contains__(self, item):
        return item in self.handle

    def __len__(self):
        return len(self.handle)

    @property
    def attrs(self):
        return self.handle.attrs

    def keys(self):
        return list(self.handle.keys())

    def items(self):
        return list(self.handle.items())

    def close(self):
        if self.handle is not None:
            self.handle.close()


class FITSFileHandler(HDF5FileHandler):
    def __init__(self, filename):
        from yt.utilities.on_demand_imports import _astropy

        if isinstance(filename, _astropy.pyfits.hdu.image._ImageBaseHDU):
            self.handle = _astropy.pyfits.HDUList(filename)
        elif isinstance(filename, _astropy.pyfits.HDUList):
            self.handle = filename
        else:
            self.handle = _astropy.pyfits.open(
                filename, memmap=True, do_not_scale_image_data=True, ignore_blank=True
            )
        self._fits_files = []

    def __del__(self):
        for f in self._fits_files:
            f.close()
        del self._fits_files
        del self.handle
        self.handle = None

    def close(self):
        self.handle.close()


def valid_netcdf_classic_signature(filename):
    signature_v1 = b"CDF\x01"
    signature_v2 = b"CDF\x02"
    try:
        with open(filename, "rb") as f:
            header = f.read(4)
            return header == signature_v1 or header == signature_v2
    except Exception:
        return False


def warn_netcdf(fn):
    # There are a few variants of the netCDF format.
    classic = valid_netcdf_classic_signature(fn)
    # NetCDF-4 Classic files are HDF5 files constrained to the Classic
    # data model used by netCDF-3.
    netcdf4_classic = valid_hdf5_signature(fn) and fn.endswith((".nc", ".nc4"))
    needs_netcdf = classic or netcdf4_classic
    from yt.utilities.on_demand_imports import _netCDF4 as netCDF4

    if needs_netcdf and isinstance(netCDF4.Dataset, NotAModule):
        raise RuntimeError(
            "This appears to be a netCDF file, but the "
            "python bindings for netCDF4 are not installed."
        )


class NetCDF4FileHandler:
    def __init__(self, filename):
        self.filename = filename

    @contextmanager
    def open_ds(self, **kwargs):
        from yt.utilities.on_demand_imports import _netCDF4 as netCDF4

        ds = netCDF4.Dataset(self.filename, mode="r", **kwargs)
        yield ds
        ds.close()

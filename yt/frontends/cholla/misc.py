from yt.utilities.on_demand_imports import _h5py


class _CachedH5Openner:
    """
    A simple context manager that helps implement the idiom where data is read
    from (or written to) one or more HDF5 and we want to wait to close the
    previous HDF5 file until it is time to open a new file. This lets us avoid
    overhead in cases where we would close and then immediately reopen the
    same file.

    By using a context manager, we're able to properly cleanup in the event
    that an exception occurs.
    """

    def __init__(self, mode: str = "r"):
        self._filename = None
        self._fh = None
        self._mode = mode

    def open_fh(self, filename: str):
        if self._filename == filename:
            return self._fh
        if self._fh is not None:
            self._fh.close()
        self._fh = _h5py.File(filename, self._mode)
        self._filename = filename
        return self._fh

    def __enter__(self):
        return self

    def __exit__(self, exc, value, tb):
        if self._fh is not None:
            self._fh.close()

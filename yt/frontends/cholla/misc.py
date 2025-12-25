import os
import typing
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

# this is a hacky workaround to get _h5py.File to work in annotations. We can probably
# address this issue more robustly by directly modifying yt.utilities.on_demand_imports
if typing.TYPE_CHECKING:
    import h5py as _h5py
else:
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

    def __init__(self, mode="r"):
        self._filename = None
        self._fh = None
        self._mode = mode

    def open_fh(self, filename):
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


@dataclass(kw_only=True, slots=True, frozen=True)
class _BlockDiskMapping:
    """Contains info for mapping blockids to locations in hdf5 files

    Notes
    -----
    At the time of writing, this is primarily meant to provide a mapping for
    field data. In the future, we may initialize a separate instance to
    provide a mapping for particle data
    """

    # ``fname_template.format(blockid=...)`` produces the file containing blockid (this
    # can properly handle cases where all blocks are stored in a single file)
    fname_template: str
    # hdf5 group containing the field data
    h5_group: str
    # maps blockid to an index that select all associated data from a field-dataset
    idx_map: Mapping[int, tuple[int | slice, ...]]


def _infer_blockid_location_arr(fname_template, global_dims, arr_shape):
    # used when hdf5 files don't have an explicit "domain" group
    blockid_location_arr = np.empty(shape=tuple(int(e) for e in arr_shape), dtype="i8")
    if blockid_location_arr.size == 1:
        # primarily intended to handle the result of older concatenation scripts (it
        # also handles the case when only a single block is used, which is okay)
        blockid_location_arr[0, 0, 0] = 0
    else:  # handle distributed cholla datasets
        local_dims, rem = np.divmod(global_dims, blockid_location_arr.shape)
        assert np.all(rem == 0) and np.all(local_dims > 0)
        for blockid in range(0, blockid_location_arr.size):
            with _h5py.File(fname_template.format(blockid=blockid), "r") as f:
                tmp, rem = np.divmod(f.attrs["offset"][:], local_dims)
            assert np.all(rem == 0)  # sanity check
            idx3D = tuple(int(e) for e in tmp)
            blockid_location_arr[idx3D] = blockid
    return blockid_location_arr


def _determine_data_layout(f: _h5py.File) -> tuple[np.ndarray, _BlockDiskMapping]:
    """Determine the data layout of the snapshot

    The premise is that the basic different data formats shouldn't
    matter outside of this function."""
    filename = f.filename

    # STEP 1: infer the template for all Cholla data-files by inspecting filename
    # ===========================================================================
    # There are 2 conventions for the names of Cholla's data-files:
    #  1. "root.h5.{blockid}" is the standard format Cholla uses when writing files
    #     storing a single snapshot. Each MPI-rank will write a separate file and
    #     replace ``{blockid}`` with MPI-rank (Modern Cholla versions without MPI
    #     replace ``{blockid}`` with ``0``)
    #  2. "root.h5": is the standard format used by Cholla's concatenation scripts
    #     (older versions of Cholla without MPI also used this format to name outputs)
    inferred_fname_template, cur_filename_suffix = _infer_fname_template(filename)

    # STEP 2: Check whether the hdf5 file has a flat structure
    # ========================================================
    # Historically, we would always store datasets directly in the root group of the
    # data file. More recent concatenation scripts store no data in groups.
    flat_structure = any(not isinstance(elem, _h5py.Group) for elem in f.values())

    # STEP 3: Extract basic domain info information from the file(s)
    # ==============================================================
    has_explicit_domain_info = "domain" in f
    if has_explicit_domain_info:
        # this branch primarily handles concatenated files made with newer logic
        blockid_location_arr = f["domain/blockid_location_arr"][...]
        field_idx_map = {
            int(blockid): (i, slice(None), slice(None), slice(None))
            for i, blockid in enumerate(f["domain/stored_blockid_list"][...])
        }
        consolidated_data = len(field_idx_map) == blockid_location_arr.size
        if not consolidated_data:
            # in the near future, we may support one of the 2 cases:
            # > if (flat_structure):
            # >     _common_idx = (slice(None), slice(None), slice(None))
            # > else:
            # >     _common_idx = (0, slice(None), slice(None), slice(None))
            # > field_idx_map = defaultdict(lambda arg=_common_idx: arg)
            raise ValueError(
                "no support for reading Cholla datasets where data is distributed "
                "among files that explicitly encode domain info."
            )
    else:  # (not has_explicit_domain_info)
        # this branch covers distributed datasets (directly written by Cholla) and
        # older concatenated files.
        #
        # historically, when the dataset is concatenated (in post-processing),
        # the "nprocs" hdf5 attribute has been dropped
        blockid_location_arr = _infer_blockid_location_arr(
            fname_template=inferred_fname_template,
            global_dims=f.attrs["dims"].astype("=i8"),
            arr_shape=f.attrs.get("nprocs", np.array([1, 1, 1])).astype("=i8"),
        )
        consolidated_data = blockid_location_arr.size == 1

        def _get_common_idx():
            return (slice(None), slice(None), slice(None))

        field_idx_map = defaultdict(_get_common_idx)

    # STEP 4: Finalize the fname template
    # ===================================
    if consolidated_data:
        fname_template = filename
    elif cur_filename_suffix != 0:
        raise ValueError(  # mostly just a sanity check!
            "filename passed to yt.load for a distributed cholla dataset must "
            "end in '.0'"
        )
    else:
        fname_template = inferred_fname_template

    mapping = _BlockDiskMapping(
        fname_template=fname_template,
        h5_group="./" if flat_structure else "field",
        idx_map=field_idx_map,
    )
    return blockid_location_arr, mapping


def _infer_fname_template(filename: str) -> tuple[str, int | None]:
    """Infers the template for all Cholla data-files based on the filename
    passed to ``yt.load``.

    string from the process-id suffix, and returns both parts in a 2-tuple.

    There are 2 conventions for the names of Cholla's data-files:
      1. "root.h5.{blockid}" is the standard format Cholla uses when writing
         files storing a single snapshot. Each MPI-rank will write a separate
         file and replace ``{blockid}`` with MPI-rank (Modern Cholla versions
         without MPI replace ``{blockid}`` with ``0``)
      2. "root.h5": is the standard format used by Cholla's concatenation
         scripts (older versions of Cholla without MPI also used this format
         to name outputs)

    Returns
    -------
    template: str
        The path to the file containing a blockid is given by calling
        ``template.format(blockid=<blockid>)``. (This will work whether
        all blocks are stored in 1 file or blocks are distributed across
        files)
    cur_blockid_suffix: int or None
        The blockid specified in the suffix of ``filename``. If there isn't a
        suffix, then this will be None.
    """

    # at this time, we expect the suffix to be the minimum number of characters
    # that are necessary to represent the process id. For flexibility, we will
    # allow extra zero-padding

    _dir, _base = os.path.split(filename)
    match _base.rpartition("."):
        case ("", ".", _):  # Cholla never writes a file like this
            raise ValueError(
                f"1st character in {filename!r} is the only '.' in the file's name"
            )
        case (prefix, ".", suffix) if suffix.isdecimal():
            return os.path.join(_dir, f"{prefix}.{{blockid}}"), int(suffix)
        case _:
            return (filename, None)

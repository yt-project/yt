"""
Tests the Cholla-frontend by generating synthetic datasets

In more detail, Cholla has had a couple of historical data formats.
- see the ``ChollaDataFmt`` enumeration for a description of each format
- we choose to use synthetic datasets in case there is a future where the
  data-format changes again and we want to maintain backwards compatability
  without uploading more and more sample Cholla datasets
- the logic for creating synthetic datasets is adapted from similar logic
  used for testing the ``cholla_utils`` python package
  - the ``cholla_utils`` python package is developed within the Cholla
    repository that provides a light-weight (compared to yt) interface for
    loading datasets
  - ideally, we will try to keep the testing logic relatively consistent
    between the 2 packages
"""

import enum
import typing
from collections.abc import Sequence

import numpy as np
import pytest

import yt
from yt.frontends.cholla.misc import _CachedH5Openner
from yt.testing import requires_module

# this is a hacky workaround to get h5py.File to work in annotations. We can probably
# address this issue more robustly by directly modifying yt.utilities.on_demand_imports
if typing.TYPE_CHECKING:
    import h5py
else:
    from yt.utilities.on_demand_imports import _h5py as h5py


class ChollaDataFmt(enum.Enum):
    """Describes the format of the grid data"""

    # the format directly written by Cholla (each block is written to a separate file)
    DISTRIBUTED = (enum.auto(), False)
    # Cholla's older concatenation scripts (that are no longer available), would
    # combine all blocks into 1 giant block. The resulting generally appears as if
    # Cholla was run with a single process that evolved a single giant block of data
    LEGACY_CONCAT = (enum.auto(), True)
    # Cholla's newer concatenation scripts combine all of the data into a single file,
    # but retains the original block structure
    CONCAT = (enum.auto(), True)

    def __new__(cls, value: typing.Any, is_single_file: bool):
        # based on example from docs
        if isinstance(value, enum.auto):
            value = len(cls.__members__) + 1

        obj = object.__new__(cls)
        obj._value_ = value
        obj.is_single_file = is_single_file
        return obj

    def __repr__(self):
        # based on example from docs (when we want to hide the underlying value)
        return f"<{self.__class__.__name__}, {self.name}>"


def _generate_array(shape: tuple[int, ...], *, start: int = 0):
    # used to generate an array of unique values of a given shape
    size = np.prod(shape)
    return np.arange(start, start + size).reshape(shape).astype("f8")


def _add_standard_header_attrs(f: h5py.File):
    # we could customize this quite a bit... (but, that seems unnecessary)

    # fields that must be handled separately:
    # - particle & field headers always have:
    #   * "bounds", "dx", "domain" (handled by _write_domain_prop_attrs)
    #   * "dims", "n_fields", "nprocs"
    #   * sometimes: "offset", "dims_local" (depends on the style)
    # - particle headers may also have:
    #   * "dt_particles" (this can be different from "dt")
    #   * "t_particles" (as far as I can tell, this is the same as "t")
    #   * sometimes: "n_particles_local" (depends on the file-style)

    f.attrs["Git Commit Hash"] = np.array(["<garbage>"], dtype=object)
    f.attrs["Macro Flags"] = np.array(["<garbage>"], dtype=object)
    f.attrs["cholla"] = np.array([""], dtype=object)
    f.attrs["density_unit"] = np.array([6.76810999e-32], dtype="f8")
    f.attrs["energy_unit"] = np.array([6.47112563e-10], dtype="f8")
    f.attrs["gamma"] = np.array([1.66666667], dtype="f8")
    f.attrs["length_unit"] = np.array([3.08567758e21], dtype="f8")
    f.attrs["mass_unit"] = np.array([1.98847e33], dtype="f8")
    f.attrs["time_unit"] = np.array([3.15569e10], dtype="f8")
    f.attrs["velocity_unit"] = np.array([9.77813911e10], dtype="f8")
    f.attrs["n_step"] = np.array([0], dtype="i4")
    f.attrs["t"] = np.array([0.0], dtype="f8")
    f.attrs["dt"] = np.array([0.0], dtype="f8")


def _write_domain_prop_attrs(f: h5py.File, global_shape: tuple[int, ...]):
    # we could customize this quite a bit... (but, that seems unnecessary)

    dx = np.array([1.0 for _ in global_shape])
    domain = dx * global_shape
    f.attrs["dx"] = dx
    f.attrs["domain"] = domain
    f.attrs["bounds"] = -0.25 * domain


def _generate_files(
    root_path: str,
    nprocs: Sequence[int],
    global_shape: Sequence[int],
    *,
    data_format: ChollaDataFmt | None = None,
    field_names: Sequence[str] | None = None,
) -> tuple[dict[str, np.ndarray], str]:
    """
    Generates file(s) that emulate a dataset holding results from a
    hypothetical Cholla simulation

    Parameters
    ----------
    root_path
        Prefix of the path where the dataset is written
    nprocs
        Specifies the number of processes that the hypothetical Cholla
        simulation used
    global_shape
        Specifies the global shape of each field in the hypothetical
        Cholla simulation (each MPI process would be responsible for evolving
        a subsection of the global shape).
    data_format
        Specifies the format of the dataset
    field_names
        The names of the fields to include in the output dataset

    Returns
    -------
    global_arrays: dict[str, np.ndarray]
        A dictionary of the global concatenated fields
    root_fname: str
        Path to one of the files. For distributed datasets this is
        always process 0.
    """
    # check and sanitize arguments

    if isinstance(field_names, str):
        raise TypeError("field_names can't be a string")
    elif field_names is None:
        field_names = ["density"]
    elif len(field_names) != len(set(field_names)):
        raise ValueError("field_names must hold unique names")

    if len(nprocs) != 3:
        raise ValueError("nprocs must be a 3 element array")
    elif any(int(e) != e for e in nprocs):
        raise ValueError("nprocs must contain integers")
    elif any(e < 1 for e in nprocs):
        raise ValueError("nprocs must contain positive values")
    else:
        nprocs = tuple(int(e) for e in nprocs)

    if len(global_shape) != 3:
        raise ValueError("global_shape must be a 3 element array")
    elif any(int(e) != e for e in global_shape):
        raise ValueError("global_shape must contain integers")
    elif any(e < 1 for e in global_shape):
        raise ValueError("global_shape must contain positive values")
    else:
        global_shape = tuple(int(e) for e in global_shape)

    # infer the shape of each block (and perform a sanity check)
    cc_block_shape, remainder = np.divmod(global_shape, nprocs)
    if (cc_block_shape == 0).any():
        raise ValueError(
            "nprocs contains a value exceeding the corrsponding length in global_shape"
        )
    elif (remainder != 0).any():
        raise ValueError(
            "at least 1 element of global_shape isn't evenly divisible by nprocs"
        )

    if data_format is ChollaDataFmt.LEGACY_CONCAT:
        cc_block_shape = np.array(global_shape)
        nprocs = (1, 1, 1)

    # construct the global arrays
    global_arrays = {}
    start_offset = 0
    for name in field_names:
        arr = _generate_array(global_shape, start=start_offset)
        global_arrays[name] = arr
        start_offset = arr.max() + 1

    # prepare to creating the files
    blockid_location_arr = np.arange(np.prod(nprocs)).reshape(nprocs)
    match data_format:
        case ChollaDataFmt.DISTRIBUTED:
            fname_template = f"{root_path}.h5.{{blockid:d}}"
            field_grp = "/"
            field_dset_shape = tuple(cc_block_shape)
        case ChollaDataFmt.LEGACY_CONCAT:
            fname_template = f"{root_path}.h5"
            field_grp = "/"
            field_dset_shape = global_shape
        case ChollaDataFmt.CONCAT:
            fname_template = f"{root_path}.h5"
            field_grp = "field"
            field_dset_shape = (np.prod(nprocs),) + tuple(cc_block_shape)
        case _:
            raise RuntimeError(f"unknown data format: {data_format}")

    # actually create the file
    with _CachedH5Openner(mode="w") as h5_context_manager:
        for idx3d, blockid in np.ndenumerate(blockid_location_arr):
            f = h5_context_manager.open_fh(fname_template.format(blockid=blockid))

            # selects the region of a global arrays relevant for the current block
            src_slc = (
                slice(idx3d[0] * cc_block_shape[0], (idx3d[0] + 1) * cc_block_shape[0]),
                slice(idx3d[1] * cc_block_shape[1], (idx3d[1] + 1) * cc_block_shape[1]),
                slice(idx3d[2] * cc_block_shape[2], (idx3d[2] + 1) * cc_block_shape[2]),
            )

            # determine the region of the output dataset relevant for the current block
            # and write any extra output-specific metadata
            match data_format:
                case ChollaDataFmt.DISTRIBUTED:
                    dst_idx = (...,)

                    f.attrs["offset"] = np.array([int(slc.start) for slc in src_slc])
                    f.attrs["dims_local"] = cc_block_shape

                case ChollaDataFmt.LEGACY_CONCAT:
                    dst_idx = (...,)

                case ChollaDataFmt.CONCAT:
                    dst_idx = (blockid, ...)

                    if blockid == 0:
                        f.create_group("domain")
                        f["domain"]["blockid_location_arr"] = blockid_location_arr
                        f["domain"]["stored_blockid_list"] = np.arange(
                            blockid_location_arr.size
                        )

                        f.create_group("field")

                case _:
                    raise RuntimeError(f"unknown data format: {data_format}")

            if (blockid == 0) or (not data_format.is_single_file):
                # write some common metadata
                f.attrs["dims"] = np.array(global_shape)
                f.attrs["nprocs"] = np.array(nprocs)
                f.attrs["n_fields"] = np.array([len(field_names)])
                _add_standard_header_attrs(f)
                _write_domain_prop_attrs(f, global_shape=global_shape)
                # create the datasets that will hold the fields
                for field_name in field_names:
                    f[field_grp].create_dataset(
                        name=field_name, shape=field_dset_shape, dtype="f8"
                    )

            # actually record the field data
            for field_name in field_names:
                f[field_grp][field_name][dst_idx] = global_arrays[field_name][src_slc]
    return global_arrays, fname_template.format(blockid=0)


_CASES = [
    {
        "nprocs": (1, 1, 1),
        "global_shape": (8, 8, 8),
        "data_format": ChollaDataFmt.DISTRIBUTED,
    },
    {
        "nprocs": (2, 2, 2),
        "global_shape": (4, 16, 8),
        "data_format": ChollaDataFmt.DISTRIBUTED,
    },
    {
        "nprocs": (1, 4, 2),
        "global_shape": (4, 16, 8),
        "data_format": ChollaDataFmt.DISTRIBUTED,
    },
    # there no point going through lots of varieties of ChollaDataFmt.LEGACY_CONCAT
    # -> the files always look very similar to each other
    {
        "nprocs": (2, 2, 2),
        "global_shape": (4, 16, 8),
        "data_format": ChollaDataFmt.LEGACY_CONCAT,
    },
    # it's definitely worth checking ChollaDataFmt.LEGACY_CONCAT when there is only
    # 1 process as well as when there are multiple processes
    {
        "nprocs": (1, 1, 1),
        "global_shape": (8, 8, 8),
        "data_format": ChollaDataFmt.CONCAT,
    },
    {
        "nprocs": (2, 2, 2),
        "global_shape": (4, 16, 8),
        "data_format": ChollaDataFmt.CONCAT,
    },
]


@requires_module("h5py")
@pytest.mark.parametrize("kwargs", _CASES)
def test_load(tmp_path, kwargs):
    # generate a synthetic dataset and make sure that the loaded values are correct

    # Step 1: create the synthetic dataset
    # -> global_arr is a dict that maps each field_name to an array that holding the
    #    fully concatenated array that holds the expected field values
    # -> root_fname is the path that should be passed to yt.load
    global_arr, root_fname = _generate_files(
        root_path=str(tmp_path / "0"), field_names=["density", "momentum_x"], **kwargs
    )

    # Step 2: load the dataset and build a covering grid spanning the full domain
    # -> accessing a field from the covering grid should return a field with the same
    #    shape as accessing the corresponding field in global_arr
    ds = yt.load(root_fname)
    grid = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )

    # Step 3: actually compare the loaded field values against our expectations
    for field_name, expected_arr in global_arr.items():
        np.testing.assert_equal(
            actual=grid["cholla", field_name].ndview,
            desired=expected_arr,
            err_msg=f"there was an issue with loading the {field_name!r} field",
        )

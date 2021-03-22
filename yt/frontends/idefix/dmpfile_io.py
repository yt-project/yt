import re
import struct
from typing import List

import numpy as np

# hardcoded in idefix
HEADERSIZE = 128
NAMESIZE = 16

SIZE_CHAR = 1
SIZE_INT = 4
# emulating C++
# enum DataType {DoubleType, SingleType, IntegerType};
DTYPES = ["d", "f", "i"]


def read_null_terminated_string(fh, maxsize=NAMESIZE):
    """Read maxsize * SIZE_CHAR bytes, but only parse non-null characters."""
    b = fh.read(maxsize * SIZE_CHAR)
    s = b.decode("utf-8", errors="backslashreplace")
    s = s.split("\x00", maxsplit=1)[0]
    return s


def read_next_field_properties(fh):
    """Emulate Idefix's OutputDump::ReadNextFieldProperty"""
    field_name = read_null_terminated_string(fh)

    fmt = "=i"
    dtype = DTYPES[struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]]
    ndim = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    if ndim > 3:
        raise ValueError(ndim)
    fmt = f"={ndim}i"
    dim = np.array(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    return field_name, dtype, ndim, dim


def read_chunk(fh, ndim: int, dim: List[int], dtype, is_scalar=False, skip_data=False):
    assert ndim == len(dim)
    fmt = f"={np.product(dim)}{dtype}"
    size = struct.calcsize(fmt)
    if skip_data:
        fh.seek(size, 1)
        return None
    data = struct.unpack(fmt, fh.read(size))

    # note: this reversal may not be desirable in general
    if is_scalar:
        return data[0]

    data = np.reshape(data, dim, order="F")
    return data


def read_serial(fh, ndim: int, dim: List[int], dtype, is_scalar=False):
    """Emulate Idefix's OutputDump::ReadSerial"""
    assert ndim == 1  # corresponds to an error raised in IDEFIX
    return read_chunk(fh, ndim=ndim, dim=dim, dtype=dtype, is_scalar=is_scalar)


def read_distributed(fh, dim, skip_data=False):
    """Emulate Idefix's OutputDump::ReadDistributed"""
    # note: OutputDump::ReadDistributed only read doubles
    # because chucks written in integers are small enough
    # that parallelization is counter productive.
    # This a design choice on idefix's size.
    return read_chunk(fh, ndim=len(dim), dim=dim, dtype="d", skip_data=skip_data)


# The following functions are originally designed for yt


def read_header(filename):
    with open(filename, "rb") as fh:
        header = read_null_terminated_string(fh, maxsize=HEADERSIZE)
    return header


def parse_fields_index(fh):
    """
    Go over a dumpfile, parse bytes offsets associated with each field.
    """
    field_index = {}

    # skip header
    fh.seek(HEADERSIZE)
    # skip grid properties
    for _ in range(9):
        _field_name, dtype, ndim, dim = read_next_field_properties(fh)
        read_serial(fh, ndim, dim, dtype)

    while True:
        offset = fh.tell()
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        if not re.match("^V[cs]-", field_name):
            break
        field_index[field_name] = offset
        read_distributed(fh, dim, skip_data=True)

    return field_index


def read_single_field(fh, field_offset):
    fh.seek(field_offset)
    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    data = read_distributed(fh, dim)
    return data


def read_idefix_dmpfile(filename, skip_data=False):
    with open(filename, "rb") as fh:
        return read_idefix_dump_from_buffer(fh, skip_data)


def read_idefix_dump_from_buffer(fh, skip_data=False):

    # skip header
    fh.seek(HEADERSIZE)

    fprops = {}
    fdata = {}
    for _ in range(9):
        # read grid properties
        # (cell centers, left and right edges in 3D -> 9 arrays)
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        data = read_serial(fh, ndim, dim, dtype)
        fprops[field_name] = dtype, ndim, dim
        fdata[field_name] = data

    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    while field_name != "eof":
        # note that this could likely be implemented using a call to
        # `iter` with a sentinel value, to the condition that read_next_field_properties
        # would be splitted into 2 parts (I don't the sentinel pattern works with tuples)
        fprops[field_name] = dtype, ndim, dim
        if field_name.startswith("Vc-") or field_name.startswith("Vs-"):
            data = read_distributed(fh, dim, skip_data=skip_data)
        else:
            is_scalar = ndim == 1 and dim[0] == 1
            is_scalar &= field_name not in ("x1", "x2", "x3")
            data = read_serial(fh, ndim, dim, dtype, is_scalar=is_scalar)
        fdata[field_name] = data
        field_name, dtype, ndim, dim = read_next_field_properties(fh)

    return fprops, fdata

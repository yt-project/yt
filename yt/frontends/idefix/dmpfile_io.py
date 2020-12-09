import struct
from typing import List

import numpy as np

# hardcoded in idefix
HEADERSIZE = 128
NAMESIZE = 16

SIZE_CHAR = 1
SIZE_INT = 4
# emulating CPP
# enum DataType {DoubleType, SingleType, IntegerType};
DTYPES = ["d", "f", "i"]

## the following methods are translations from c++ to Python


def read_str(fh, size=NAMESIZE):
    fmt = "=" + size * "c"
    raw_cstring64 = iter(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    c = next(raw_cstring64)
    s = ""
    while r"\x" not in c.__str__():  # todo: better condition here
        # emulating (poorly) std::strlen
        # read NAMESIZE * SIZE_CHAR bytes, but only parse non-null characters
        s += c.decode()
        c = next(raw_cstring64)
    return s


def read_next_field_properties(fh):

    field_name = read_str(fh)

    fmt = "=i"
    dtype = DTYPES[struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]]
    ndim = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    if ndim > 3:
        raise ValueError(ndim)
    fmt = "=" + ndim * "i"
    dim = np.array(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    return field_name, dtype, ndim, dim


def read_chunk(fh, ndim: int, dim: List[int], dtype, is_scalar=False, skip_data=False):
    assert ndim == len(dim)
    fmt = "=" + np.product(dim) * dtype
    size = struct.calcsize(fmt)
    if skip_data:
        fh.seek(size, 1)
        return
    data = struct.unpack(fmt, fh.read(size))

    # note: this reversal may not be desirable in general
    if is_scalar:
        return data[0]

    data = np.reshape(data, dim[::-1])
    return data


def read_serial(fh, ndim: int, dim: List[int], dtype, is_scalar=False):
    assert ndim == 1  # corresponds to an error raised in IDEFIX
    return read_chunk(fh, ndim=ndim, dim=dim, dtype=dtype, is_scalar=is_scalar)


def read_distributed(fh, dim, skip_data=False):
    # note: OutputDump::ReadDistributed only read doubles
    # because chucks written in integers are small enough
    # that parallelization is counter productive.
    # This a design choice on idefix's size.
    return read_chunk(fh, ndim=len(dim), dim=dim, dtype="d", skip_data=skip_data)


def read_header(filepath_or_buffer):
    if "read" in filepath_or_buffer.__dir__():
        fh = filepath_or_buffer
        closeme = False
    else:
        fh = open(filepath_or_buffer, "rb")
        closeme = True
    fh.seek(0)
    header = read_str(fh, size=HEADERSIZE)
    if closeme:
        fh.close()
    return header


def read_idefix_dmpfile(filepath_or_buffer, skip_data=False):
    fprops = {}
    fdata = {}
    if "read" in filepath_or_buffer.__dir__():
        fh = filepath_or_buffer
        closeme = False
    else:
        fh = open(filepath_or_buffer, "rb")
        closeme = True

    fh.seek(0)

    # skip header
    read_header(fh)

    for _ in range(9):
        # read grid properties
        # (cell centers, left and right edges in 3D -> 9 arrays)
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        data = read_serial(fh, ndim, dim, dtype)
        fprops.update({field_name: (dtype, ndim, dim)})
        fdata.update({field_name: data})

    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    while field_name != "eof":
        fprops.update({field_name: (dtype, ndim, dim)})
        if field_name.startswith("Vc-") or field_name.startswith("Vs-"):
            data = read_distributed(fh, dim, skip_data=skip_data)
        else:
            is_scalar = ndim == 1 and dim[0] == 1
            is_scalar &= field_name not in ("x1", "x2", "x3")
            data = read_serial(fh, ndim, dim, dtype, is_scalar=is_scalar)
        fdata.update({field_name: data})
        field_name, dtype, ndim, dim = read_next_field_properties(fh)

    if closeme:
        fh.close()
    return fprops, fdata

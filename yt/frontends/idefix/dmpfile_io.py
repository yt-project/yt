import struct
from typing import List

import numpy as np

NAMESIZE = 16  # maximum size of the name array (hardcoded in idefix)
SIZE_CHAR = 1
SIZE_INT = 4
# emulating CPP
# enum DataType {DoubleType, SingleType, IntegerType};
DTYPES = [np.float64, np.float32, np.int32]

## the following methods are translations from c++ to Python


def read_next_field_properties(fh):
    fmt = "=" + NAMESIZE * "c"
    raw_cstring64 = iter(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    c = next(raw_cstring64)
    field_name = ""
    while r"\x" not in c.__str__():  # todo: better condition here
        # emulating (poorly) std::strlen
        # read NAMESIZE * SIZE_CHAR bytes, but only parse non-null characters
        field_name += c.decode()
        c = next(raw_cstring64)

    fmt = "=i"
    dtype = DTYPES[struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]]
    ndim = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    if ndim > 3:
        raise ValueError(ndim)
    fmt = "=" + ndim * "i"
    dim = np.array(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    return field_name, dtype, ndim, dim


def read_serial(fh, ndim: int, dim: List[int], dtype):
    assert ndim == 1  # corresponds to an error raised in IDEFIX
    # this will need revision if this assertion ever becomes obsolete
    # note that `dim` is in general an array with size `ndim`
    # for now I'll assume `dim` is a single int
    fmt = (
        "=" + np.product(dim) * {np.float64: "d", np.float32: "f", np.int32: "i"}[dtype]
    )
    data = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    return data


def read_distributed(fh, dim):
    # note: OutputDump::ReadDistributed only read doubles
    # because chucks written in integers are small enough
    # that parallelization is counter productive.
    # This a design choice on idefix's size.

    fmt = "=" + np.product(dim) * "d"
    data = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))

    # note: this reversal may not be desirable in general
    data = np.reshape(data, dim[::-1])
    return data


def read_idefix_dmpfile(filepath_or_buffer):
    fprops = {}
    fdata = {}
    if "read" in filepath_or_buffer.__dir__():
        fh = filepath_or_buffer
        closeme = False
    else:
        fh = open(filepath_or_buffer, "rb")
        closeme = True

    for _ in range(3):
        # read grid properties
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        data = read_serial(fh, ndim, dim, dtype)
        fprops.update({field_name: (dtype, ndim, dim)})
        fdata.update({field_name: data})

    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    while field_name != "eof":
        fprops.update({field_name: (dtype, ndim, dim)})
        if field_name.startswith("Vc-") or field_name.startswith("Vs-"):
            data = read_distributed(fh, dim)
        else:
            data = read_serial(fh, ndim, dim, dtype)
        fdata.update({field_name: data})
        field_name, dtype, ndim, dim = read_next_field_properties(fh)

    if closeme:
        fh.close()
    return fprops, fdata

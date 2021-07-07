import io
import os
import struct

import numpy as np


def read_attrs(f, attrs, endian="="):
    r"""This function accepts a file pointer and reads from that file pointer
    according to a definition of attributes, returning a dictionary.

    Fortran unformatted files provide total bytesize at the beginning and end
    of a record.  By correlating the components of that record with attribute
    names, we construct a dictionary that gets returned.  Note that this
    function is used for reading sequentially-written records.  If you have
    many written that were written simultaneously, see read_record.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    attrs : iterable of iterables
        This object should be an iterable of one of the formats:
        [ (attr_name, count, struct type), ... ].
        [ ((name1,name2,name3),count, vector type]
        [ ((name1,name2,name3),count, 'type type type']
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    values : dict
        This will return a dict of iterables of the components of the values in
        the file.

    Examples
    --------

    >>> header = [("ncpu", 1, "i"), ("nfiles", 2, "i")]
    >>> f = open("fort.3", "rb")
    >>> rv = read_attrs(f, header)
    """
    vv = {}
    net_format = endian
    for _a, n, t in attrs:
        for end in "@=<>":
            t = t.replace(end, "")
        net_format += "".join(["I"] + ([t] * n) + ["I"])
    size = struct.calcsize(net_format)
    vals = list(struct.unpack(net_format, f.read(size)))
    vv = {}
    for a, n, t in attrs:
        for end in "@=<>":
            t = t.replace(end, "")
        if isinstance(a, tuple):
            n = len(a)
        s1 = vals.pop(0)
        v = [vals.pop(0) for i in range(n)]
        s2 = vals.pop(0)
        if s1 != s2:
            size = struct.calcsize(endian + "I" + "".join(n * [t]) + "I")
            raise OSError(
                "An error occurred while reading a Fortran record. "
                "Got a different size at the beginning and at the "
                "end of the record: %s %s",
                s1,
                s2,
            )
        if n == 1:
            v = v[0]
        if isinstance(a, tuple):
            if len(a) != len(v):
                raise OSError(
                    "An error occurred while reading a Fortran "
                    "record. Record length is not equal to expected "
                    "length: %s %s",
                    len(a),
                    len(v),
                )
            for k, val in zip(a, v):
                vv[k] = val
        else:
            vv[a] = v
    return vv


def read_cattrs(f, attrs, endian="="):
    r"""This function accepts a file pointer to a C-binary file and reads from
    that file pointer according to a definition of attributes, returning a
    dictionary.

    This function performs very similarly to read_attrs, except it does not add
    on any record padding.  It is thus useful for using the same header types
    as in read_attrs, but for C files rather than Fortran.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    attrs : iterable of iterables
        This object should be an iterable of one of the formats:
        [ (attr_name, count, struct type), ... ].
        [ ((name1,name2,name3),count, vector type]
        [ ((name1,name2,name3),count, 'type type type']
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    values : dict
        This will return a dict of iterables of the components of the values in
        the file.

    Examples
    --------

    >>> header = [("ncpu", 1, "i"), ("nfiles", 2, "i")]
    >>> f = open("cdata.bin", "rb")
    >>> rv = read_cattrs(f, header)
    """
    vv = {}
    net_format = endian
    for _a, n, t in attrs:
        for end in "@=<>":
            t = t.replace(end, "")
        net_format += "".join([t] * n)
    size = struct.calcsize(net_format)
    vals = list(struct.unpack(net_format, f.read(size)))
    vv = {}
    for a, n, t in attrs:
        for end in "@=<>":
            t = t.replace(end, "")
        if isinstance(a, tuple):
            n = len(a)
        v = [vals.pop(0) for i in range(n)]
        if n == 1:
            v = v[0]
        if isinstance(a, tuple):
            if len(a) != len(v):
                raise OSError(
                    "An error occurred while reading a Fortran "
                    "record. Record length is not equal to expected "
                    "length: %s %s",
                    len(a),
                    len(v),
                )

            for k, val in zip(a, v):
                vv[k] = val
        else:
            vv[a] = v
    return vv


def read_vector(f, d, endian="="):
    r"""This function accepts a file pointer and reads from that file pointer
    a vector of values.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    d : data type
        This is the datatype (from the struct module) that we should read.
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    tr : numpy.ndarray
        This is the vector of values read from the file.

    Examples
    --------

    >>> f = open("fort.3", "rb")
    >>> rv = read_vector(f, "d")
    """
    pad_fmt = f"{endian}I"
    pad_size = struct.calcsize(pad_fmt)
    vec_len = struct.unpack(pad_fmt, f.read(pad_size))[0]  # bytes
    vec_fmt = f"{endian}{d}"
    vec_size = struct.calcsize(vec_fmt)
    if vec_len % vec_size != 0:
        raise OSError(
            "An error occurred while reading a Fortran record. "
            "Vector length is not compatible with data type: %s %s",
            vec_len,
            vec_size,
        )
    vec_num = int(vec_len / vec_size)
    if isinstance(f, io.IOBase):
        tr = np.frombuffer(f.read(vec_len), vec_fmt, count=vec_num)
    else:
        tr = np.frombuffer(f, vec_fmt, count=vec_num)
    vec_len2 = struct.unpack(pad_fmt, f.read(pad_size))[0]
    if vec_len != vec_len2:
        raise OSError(
            "An error occurred while reading a Fortran record. "
            "Got a different size at the beginning and at the "
            "end of the record: %s %s",
            vec_len,
            vec_len2,
        )
    return tr


def skip(f, n=1, endian="="):
    r"""This function accepts a file pointer and skips a Fortran unformatted
    record. Optionally check that the skip was done correctly by checking
    the pad bytes.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    n : int
        Number of records to skip.
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    skipped: The number of elements in the skipped array

    Examples
    --------

    >>> f = open("fort.3", "rb")
    >>> skip(f, 3)
    """
    skipped = np.zeros(n, dtype=np.int32)
    fmt = endian + "I"
    fmt_size = struct.calcsize(fmt)
    for i in range(n):
        size = f.read(fmt_size)
        s1 = struct.unpack(fmt, size)[0]
        f.seek(s1 + fmt_size, os.SEEK_CUR)
        s2 = struct.unpack(fmt, size)[0]
        if s1 != s2:
            raise OSError(
                "An error occurred while reading a Fortran record. "
                "Got a different size at the beginning and at the "
                "end of the record: %s %s",
                s1,
                s2,
            )

        skipped[i] = s1 / fmt_size
    return skipped


def peek_record_size(f, endian="="):
    r"""This function accept the file handle and returns
    the size of the next record and then rewinds the file
    to the previous position.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    Number of bytes in the next record
    """
    pos = f.tell()
    s = struct.unpack(">i", f.read(struct.calcsize(">i")))
    f.seek(pos)
    return s[0]


def read_record(f, rspec, endian="="):
    r"""This function accepts a file pointer and reads from that file pointer
    a single "record" with different components.

    Fortran unformatted files provide total bytesize at the beginning and end
    of a record.  By correlating the components of that record with attribute
    names, we construct a dictionary that gets returned.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    rspec : iterable of iterables
        This object should be an iterable of the format [ (attr_name, count,
        struct type), ... ].
    endian : str
        '=' is native, '>' is big, '<' is little endian

    Returns
    -------
    values : dict
        This will return a dict of iterables of the components of the values in
        the file.

    Examples
    --------

    >>> header = [("ncpu", 1, "i"), ("nfiles", 2, "i")]
    >>> f = open("fort.3", "rb")
    >>> rv = read_record(f, header)
    """
    vv = {}
    net_format = endian + "I"
    for _a, n, t in rspec:
        t = t if len(t) == 1 else t[-1]
        net_format += f"{n}{t}"
    net_format += "I"
    size = struct.calcsize(net_format)
    vals = list(struct.unpack(net_format, f.read(size)))
    s1, s2 = vals.pop(0), vals.pop(-1)
    if s1 != s2:
        raise OSError(
            "An error occurred while reading a Fortran record. Got "
            "a different size at the beginning and at the end of "
            "the record: %s %s",
            s1,
            s2,
        )
    pos = 0
    for a, n, _t in rspec:
        vv[a] = vals[pos : pos + n]
        pos += n
    return vv

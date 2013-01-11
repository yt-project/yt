"""
Utilities for reading Fortran files.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import struct
import numpy as np
import os

def read_attrs(f, attrs):
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
        This object should be an iterable of the format [ (attr_name, count,
        struct type), ... ].

    Returns
    -------
    values : dict
        This will return a dict of iterables of the components of the values in
        the file.

    Examples
    --------

    >>> header = [ ("ncpu", 1, "i"), ("nfiles", 2, "i") ]
    >>> f = open("fort.3", "rb")
    >>> rv = read_attrs(f, header)
    """
    vv = {}
    net_format = "="
    for a, n, t in attrs:
        net_format += "".join(["I"] + ([t] * n) + ["I"])
    size = struct.calcsize(net_format)
    vals = list(struct.unpack(net_format, f.read(size)))
    vv = {}
    for a, b, n in attrs:
        s1 = vals.pop(0)
        v = [vals.pop(0) for i in range(b)]
        s2 = vals.pop(0)
        if s1 != s2:
            size = struct.calcsize("=I" + "".join(b*[n]) + "I")
            print "S1 = %s ; S2 = %s ; %s %s %s = %s" % (
                    s1, s2, a, b, n, size)
            raise RuntimeError
        assert(s1 == s2)
        if b == 1: v = v[0]
        vv[a] = v
    return vv

def read_vector(f, d):
    r"""This function accepts a file pointer and reads from that file pointer
    a vector of values.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    d : data type
        This is the datatype (from the struct module) that we should read.

    Returns
    -------
    tr : numpy.ndarray
        This is the vector of values read from the file.

    Examples
    --------

    >>> f = open("fort.3", "rb")
    >>> rv = read_vector(f, 'd')
    """
    fmt = "=I"
    ss = struct.unpack(fmt, f.read(struct.calcsize(fmt)))[0]
    ds = struct.calcsize("=%s" % d)
    if ss % ds != 0:
        print "fmt = '%s' ; ss = %s ; ds = %s" % (fmt, ss, ds)
        raise RuntimeError
    count = ss / ds
    tr = np.fromstring(f.read(np.dtype(d).itemsize*count), d, count)
    vec = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    assert(vec[-1] == ss)
    return tr

def skip(f, n = 1):
    r"""This function accepts a file pointer and skips a Fortran unformatted
    record.

    Parameters
    ----------
    f : File object
        An open file object.  Should have been opened in mode rb.
    n : int
        Number of records to skip.

    Examples
    --------

    >>> f = open("fort.3", "rb")
    >>> skip(f, 3)
    """
    for i in range(n):
        fmt = "=I"
        ss = struct.unpack(fmt, f.read(struct.calcsize(fmt)))[0]
        f.seek(ss + struct.calcsize("=I"), os.SEEK_CUR)

def read_record(f, rspec):
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

    Returns
    -------
    values : dict
        This will return a dict of iterables of the components of the values in
        the file.

    Examples
    --------

    >>> header = [ ("ncpu", 1, "i"), ("nfiles", 2, "i") ]
    >>> f = open("fort.3", "rb")
    >>> rv = read_record(f, header)
    """
    vv = {}
    net_format = "=I" + "".join(["%s%s" % (n, t) for a, n, t in rspec]) + "I"
    size = struct.calcsize(net_format)
    vals = list(struct.unpack(net_format, f.read(size)))
    vvv = vals[:]
    s1, s2 = vals.pop(0), vals.pop(-1)
    if s1 != s2:
        print "S1 = %s ; S2 = %s ; SIZE = %s"
        raise RuntimeError
    pos = 0
    for a, n, t in rspec:
        vv[a] = vals[pos:pos+n]
        pos += n
    return vv, vvv


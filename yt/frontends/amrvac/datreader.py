# THOSE FUNCTIONS COME FROM $AMRVAC_DIR/tools/python/dat_reader.py
# ================================================================

# Tools to read in MPI-AMRVAC 2.0 snapshots (.dat files). Works for 1d, 2d and
# 3d data. For a usage example see the bottom of the file.
#
# Author: Jannis Teunissen
# This file should work in both Python 2 and 3

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import argparse
import struct
import numpy as np
import sys

# Size of basic types (in bytes)
size_logical = 4
size_int = 4
size_double = 8
name_len = 16

# For un-aligned data, use '=' (for aligned data set to '')
align = '='

def get_header(dat):
    """Read header from an MPI-AMRVAC 2.0 snapshot. Argument 'dat' should be a file
    opened in binary mode.
    """

    dat.seek(0)
    h = {}

    fmt = align + 'i'
    [h['version']] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    if h['version'] < 3:
        print("Unsupported .dat file version: ", h['version'])
        sys.exit()

    # Read scalar data at beginning of file
    fmt = align + 9 * 'i' + 'd'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    # Read min/max coordinates
    fmt = align + h['ndim'] * 'd'
    h['xmin'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
    h['xmax'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read domain and block size (in number of cells)
    fmt = align + h['ndim'] * 'i'
    h['domain_nx'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
    h['block_nx'] = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = align + name_len * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w_names.append(b''.join(hdr).strip().decode())
    h['w_names'] = w_names

    # Read physics type
    fmt = align + name_len * 'c'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['physics_type'] = b''.join(hdr).strip().decode()

    # Read number of physics-defined parameters
    fmt = align + 'i'
    [n_pars] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # First physics-parameter values are given, then their names
    fmt = align + n_pars * 'd'
    vals = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    fmt = align + n_pars * name_len * 'c'
    names = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [b''.join(names[i:i+name_len]).strip().decode()
             for i in range(0, len(names), name_len)]

    # Store the values corresponding to the names
    for val, name in zip(vals, names):
        h[name] = val
    return h

def get_block_data(dat):
    """Read block data from an MPI-AMRVAC 2.0 snapshot. Argument 'dat' should be a
    file opened in binary mode.
    """

    dat.seek(0)
    h = get_header(dat)
    nw = h['nw']
    block_nx = np.array(h['block_nx'])
    domain_nx = np.array(h['domain_nx'])
    xmax = np.array(h['xmax'])
    xmin = np.array(h['xmin'])
    nleafs = h['nleafs']
    nparents = h['nparents']

    # Read tree info. Skip 'leaf' array
    dat.seek(h['offset_tree'] + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = align + nleafs * 'i'
    block_lvls = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = align + nleafs * h['ndim'] * 'i'
    block_ixs = np.reshape(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))),
        [nleafs, h['ndim']])

    # Determine coarse grid spacing
    dx0 = (xmax-xmin) / domain_nx

    ### Start reading data blocks
    dat.seek(h['offset_blocks'])

    blocks = []

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        # Read number of ghost cells
        fmt = align + h['ndim'] * 'i'
        gc_lo = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))
        gc_hi = np.array(struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

        # Read actual data
        block_shape = np.append(gc_lo + block_nx + gc_hi, nw)
        fmt = align + np.prod(block_shape) * 'd'
        d = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w = np.reshape(d, block_shape, order='F') # Fortran ordering

        b = {}
        b['lvl'] = lvl
        b['ix'] = ix
        b['w'] = w
        blocks.append(b)
    return blocks


def get_uniform_data(dat):
    """Read block data from an MPI-AMRVAC 2.0 snapshot and convert to a uniform
    grid. Argument 'dat' should be a file opened in binary mode.
    """

    h = get_header(dat)
    blocks = get_block_data(dat)

    # Check if grid is uniformly refined
    refined_nx = 2**(h['levmax']-1) * h['domain_nx']
    nleafs_uniform = np.prod(refined_nx/h['block_nx'])

    if (h['nleafs'] == nleafs_uniform):
        domain_shape = np.append(refined_nx, h['nw'])
        d = np.zeros(domain_shape, order='F')

        for b in blocks:
            i0 = (b['ix'] - 1) * h['block_nx']
            i1 = i0 + h['block_nx']
            if h['ndim'] == 1:
                d[i0[0]:i1[0], :] = b['w']
            elif h['ndim'] == 2:
                d[i0[0]:i1[0], i0[1]:i1[1], :] = b['w']
            elif h['ndim'] == 3:
                d[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2], :] = b['w']
        return d
    else:
        sys.exit('Data in .dat file is not uniformly refined')
        return None


### The following is a refactoring of above functions by Clément Robert

def get_block_info(dat):
    """dat is assumed to be an open datfile

    This is part of get_block_data(), but only returns block levels and indices.
    This can be used as the "first pass" data reading required by YT's interface.
    """
    dat.seek(0)
    header = get_header(dat)
    nleafs = header['nleafs']
    nparents = header['nparents']

    # Read tree info. Skip 'leaf' array
    dat.seek(header['offset_tree'] + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = align + nleafs * 'i'
    block_lvls = np.array(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))))

    # Read block indices
    fmt = align + nleafs * header['ndim'] * 'i'
    block_ixs = np.reshape(
        struct.unpack(fmt, dat.read(struct.calcsize(fmt))),
        [nleafs, header['ndim']])

    return block_lvls, block_ixs

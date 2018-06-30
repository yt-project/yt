"""
Miscellaneous functions that are Enzo-P-specific



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import \
    ensure_tuple

def bdecode(block):
    """
    Decode a block descriptor to get its left and right sides and level.

    A block string consisting of (0, 1), with optionally one colon. The
    number of digits after the colon is the refinemenet level. The combined
    digits denote the binary representation of the left edge.
    """

    if ":" in block:
        level = len(block) - block.find(":") - 1
    else:
        level = 0
    bst   = block.replace(":", "")
    d     = float(2**len(bst))
    left  = int(bst, 2)
    right = left + 1
    left  /= d
    right /= d
    return level, left, right

def get_block_string_and_dim(block, min_dim=3):
    mybs = block[1:].split("_")
    dim = max(len(mybs), min_dim)
    return mybs, dim

def get_block_level(block):
    if ":" in block:
        l = block.find(":")
    else:
        l = len(block)
    return l

def get_block_info(block, min_dim=3):
    mybs, dim = get_block_string_and_dim(
        block, min_dim=min_dim)
    left = np.zeros(dim)
    right = np.ones(dim)
    for i, myb in enumerate(mybs):
        level, left[i], right[i] = bdecode(myb)
    return level, left, right

def get_root_blocks(block, min_dim=3):
    mybs, dim = get_block_string_and_dim(
        block, min_dim=min_dim)
    nb = np.ones(dim, dtype=int)
    for i, myb in enumerate(mybs):
        s = get_block_level(myb)
        nb[i] = 2**s
    return nb

def get_root_block_id(block, min_dim=3):
    mybs, dim = get_block_string_and_dim(
        block, min_dim=min_dim)
    rbid = np.zeros(dim, dtype=int)
    for i, myb in enumerate(mybs):
        s = get_block_level(myb)
        rbid[i] = int(myb[:s], 2)
    return rbid

def get_child_index(anc_id, desc_id):
    cid = ""
    for aind, dind in zip( anc_id.split("_"),
                          desc_id.split("_")):
        cid += dind[len(aind)]
    cid = int(cid, 2)
    return cid

def is_parent(anc_block, desc_block):
    dim = anc_block.count("_")
    if ( len(desc_block.replace(":", "")) -
         len( anc_block.replace(":", "")) ) / dim != 1:
        return False

    for aind, dind in zip( anc_block.split("_"),
                          desc_block.split("_")):
        if not dind.startswith(aind):
            return False
    return True

def nested_dict_get(pdict, keys, default=None):
    """
    Retrieve a value from a nested dict using a tuple of keys.

    If a is a dict, and a['b'] = {'c': 'd'},
    then nested_dict_get(a, ('b', 'c')) returns 'd'.
    """

    keys = ensure_tuple(keys)
    val = pdict
    for key in keys:
        if val is None:
            break
        val = val.get(key)
    if val is None:
        val = default
    return val

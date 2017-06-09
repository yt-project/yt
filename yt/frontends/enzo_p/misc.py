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

def bdecode(block):
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

def get_block_info(block, min_dim=3):
    mybs = block[1:].split("_")
    dim = max(len(mybs), min_dim)
    left = np.zeros(dim)
    right = np.ones(dim)
    for i, myb in enumerate(mybs):
        level, left[i], right[i] = bdecode(myb)
    return level, left, right

def get_root_blocks(block, min_dim=3):
    mybs = block[1:].split("_")
    dim = max(len(mybs), min_dim)
    nb = np.ones(dim, dtype=int)
    for i, myb in enumerate(mybs):
        if ":" in myb:
            s = myb.find(":")
        else:
            s = len(myb)
        nb[i] = 2**s
    return nb

def get_root_block_id(block, min_dim=3):
    mybs = block[1:].split("_")
    dim = max(len(mybs), min_dim)
    rbid = np.zeros(dim, dtype=int)
    for i, myb in enumerate(mybs):
        if ":" in myb:
            s = myb.find(":")
        else:
            s = len(myb)
        rbid[i] = int(myb[:s], 2)
    return rbid

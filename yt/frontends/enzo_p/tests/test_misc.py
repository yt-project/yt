"""
Enzo-P misc tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.frontends.enzo_p.misc import \
    get_block_info, \
    get_root_blocks, \
    get_root_block_id, \
    nested_dict_get

def get_random_block_string(max_n=64, random_state=None, level=None):
    if random_state is None:
        random_state = np.random.RandomState()

    max_l = int(np.log2(max_n))
    form = "%0" + str(max_l) + "d"
    num10 = random_state.randint(0, high=max_n)
    num2 = form % int(bin(num10)[2:])

    if level is None:
        level = random_state.randint(0, high=max_l)
    if level > 0:
        my_block = "%s:%s" % (num2[:-level], num2[-level:])
    else:
        my_block = num2
    my_block = "B" + my_block

    return num10, level, my_block

def test_get_block_info():
    rs = np.random.RandomState(45047)
    max_n = 64
    for i in range(10):
        n, l, b = get_random_block_string(
            max_n=max_n, random_state=rs)
        level, left, right = get_block_info(b, min_dim=1)
        assert level == l
        assert left == float(n) / max_n
        assert right == float(n+1) / max_n

def test_root_blocks():
    rs = np.random.RandomState(45652)
    for i in range(2, 6):
        max_n = 2**i
        n1, l1, b1 = get_random_block_string(
            max_n=max_n, random_state=rs, level=0)
        n2, l2, b2 = get_random_block_string(
            max_n=32, random_state=rs, level=0)
        block = "%s:%s" % (b1, b2[1:])

        nrb = get_root_blocks(block, min_dim=1)
        assert nrb == max_n
        rbid = get_root_block_id(block, min_dim=1)
        assert rbid == n1

def test_nested_dict_get():
    rs = np.random.RandomState(47988)
    keys = []
    my_dict = None
    for i in range(5):
        k = str(rs.randint(0, high=1000000))
        if my_dict is None:
            v = str(rs.randint(0, high=1000000))
            keys.append(v)
            my_dict = {k: v}
        else:
            my_dict = {k: my_dict}
        keys.append(k)
    keys.reverse()
    assert nested_dict_get(my_dict, keys[:-1]) == keys[-1]

    my_def = "devron"
    assert nested_dict_get(
        my_dict, "system", default=my_def) == my_def

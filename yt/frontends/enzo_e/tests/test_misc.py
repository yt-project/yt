import numpy as np

from yt.frontends.enzo_e.misc import (
    get_block_info,
    get_root_block_id,
    get_root_blocks,
    is_parent,
    nested_dict_get,
)


def get_random_block_string(max_n=64, random_state=None, level=None):
    if max_n == 1:
        assert level is None or level == 0
        return 0, 0, "B"
    elif max_n < 1:
        raise ValueError("max_n must be a positive integer")

    if random_state is None:
        random_state = np.random.RandomState()

    max_l = int(np.log2(max_n))
    form = f"%0{max_l}d"
    num10 = random_state.randint(0, high=max_n)
    num2 = form % int(bin(num10)[2:])  # the slice clips the '0b' prefix

    if level is None:
        level = random_state.randint(0, high=max_l)
    if level > 0:
        my_block = f"{num2[:-level]}:{num2[-level:]}"
    else:
        my_block = num2
    my_block = "B" + my_block

    return num10, level, my_block


def flip_random_block_bit(block, rs):
    """
    Flips a bit string in one of the block descriptors in a given block name
    """
    # split block into descriptors for each dimension
    descriptors = block[1:].split("_")

    # choose which descriptor to modify
    flippable = [i for i, descr in enumerate(descriptors) if len(descr) > 0]
    if len(flippable) == 0:  # when block in ['B', 'B_', 'B__']
        raise ValueError(f"{block} has no bits that can be flipped")
    descr_index = flippable[rs.randint(0, len(flippable))]

    # split block descriptor into left and right parts
    parts = descriptors[descr_index].split(":")
    # select the part to be modified
    if len(parts) == 1:  # block is unrefined
        part_index = 0
    elif len(parts[0]) == 0:  # The root block index can't be modified
        part_index = 1
    else:
        part_index = rs.randint(0, high=2)
    modify_part = parts[part_index]

    # flip a bit in modify_part, and return the new block name with this change
    flip_index = rs.randint(0, high=len(modify_part))
    parts[part_index] = "%s%d%s" % (
        modify_part[:flip_index],
        (int(modify_part[flip_index]) + 1) % 2,
        modify_part[flip_index + 1 :],
    )
    descriptors[descr_index] = ":".join(parts)
    return "B" + "_".join(descriptors)


def test_get_block_info():
    rs = np.random.RandomState(45047)
    max_n = 64
    for _ in range(10):
        n, l, b = get_random_block_string(max_n=max_n, random_state=rs)
        level, left, right = get_block_info(b, min_dim=1)
        assert level == l
        assert left == float(n) / max_n
        assert right == float(n + 1) / max_n

    for block in ["B", "B_", "B__"]:
        level, left, right = get_block_info(block)
        assert level == 0
        assert (left == 0.0).all()
        assert (right == 1.0).all()


def test_root_blocks():
    rs = np.random.RandomState(45652)
    for i in range(6):
        max_n = 2**i
        n1, l1, b1 = get_random_block_string(max_n=max_n, random_state=rs, level=0)
        n2, l2, b2 = get_random_block_string(max_n=32, random_state=rs, level=0)
        block = f"{b1}:{b2[1:]}"

        nrb = get_root_blocks(block, min_dim=1)
        assert nrb == max_n
        rbid = get_root_block_id(block, min_dim=1)
        assert rbid == n1


def test_is_parent():
    rs = np.random.RandomState(45652)
    for dim in [1, 2, 3]:
        for i in range(6):
            max_n = 2**i

            descriptors = []
            for _ in range(dim):
                n1, l1, b1 = get_random_block_string(
                    max_n=max_n, random_state=rs, level=0
                )
                n2, l2, b2 = get_random_block_string(max_n=32, random_state=rs, level=0)
                descriptors.append(f"{b1[1:]}:{b2[1:]}")
            block = "B" + "_".join(descriptors)
            # since b2 is computed with max_n=32 in the for-loop, block always
            # has a refined great-great-grandparent
            parent = "B" + "_".join(elem[:-1] for elem in descriptors)
            grandparent = "B" + "_".join(elem[:-2] for elem in descriptors)

            assert is_parent(parent, block)
            assert is_parent(grandparent, parent)
            assert not is_parent(grandparent, block)
            assert not is_parent(block, parent)
            assert not is_parent(flip_random_block_bit(parent, rs), block)


def test_nested_dict_get():
    rs = np.random.RandomState(47988)
    keys = []
    my_dict = None
    for _ in range(5):
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
    assert nested_dict_get(my_dict, "system", default=my_def) == my_def


def test_nested_dict_get_real_none():
    my_dict = {"a": None}
    response = nested_dict_get(my_dict, "a", default="fail")
    assert response is None

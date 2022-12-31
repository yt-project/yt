import numpy as np
from more_itertools import always_iterable


def bdecode(block):
    """
    Decode a block descriptor to get its left and right sides and level.

    A block string consisting of (0, 1), with optionally one colon. The
    number of digits after the colon is the refinement level. The combined
    digits denote the binary representation of the left edge.
    """

    if ":" in block:
        level = len(block) - block.find(":") - 1
    else:
        level = 0
    bst = block.replace(":", "")
    d = float(2 ** len(bst))
    left = int(bst, 2)
    right = left + 1
    left /= d
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
    """Decode a block name to get its left and right sides and level.

    Given a block name, this function returns the locations of the block's left
    and right edges (measured as binary fractions of the domain along each
    axis) and level.

    Unrefined blocks in the root array (which can each hold an of octree) have
    a refinement level of 0 while their ancestors (used internally by Enzo-E's
    solvers - they don't actually hold meaningful data) have negative levels.
    Because identification of negative refinement levels requires knowledge of
    the root array shape (the 'root_blocks' value specified in the parameter
    file), all unrefined blocks are assumed to have a level of 0.
    """
    mybs, dim = get_block_string_and_dim(block, min_dim=min_dim)
    left = np.zeros(dim)
    right = np.ones(dim)
    level = 0
    for i, myb in enumerate(mybs):
        if myb == "":
            continue
        level, left[i], right[i] = bdecode(myb)
    return level, left, right


def get_root_blocks(block, min_dim=3):
    mybs, dim = get_block_string_and_dim(block, min_dim=min_dim)
    nb = np.ones(dim, dtype="int64")
    for i, myb in enumerate(mybs):
        if myb == "":
            continue
        s = get_block_level(myb)
        nb[i] = 2**s
    return nb


def get_root_block_id(block, min_dim=3):
    mybs, dim = get_block_string_and_dim(block, min_dim=min_dim)
    rbid = np.zeros(dim, dtype="int64")
    for i, myb in enumerate(mybs):
        if myb == "":
            continue
        s = get_block_level(myb)
        if s == 0:
            continue
        rbid[i] = int(myb[:s], 2)
    return rbid


def get_child_index(anc_id, desc_id):
    cid = ""
    for aind, dind in zip(anc_id.split("_"), desc_id.split("_")):
        cid += dind[len(aind)]
    cid = int(cid, 2)
    return cid


def is_parent(anc_block, desc_block):
    dim = anc_block.count("_") + 1
    if (len(desc_block.replace(":", "")) - len(anc_block.replace(":", ""))) / dim != 1:
        return False

    for aind, dind in zip(anc_block.split("_"), desc_block.split("_")):
        if not dind.startswith(aind):
            return False
    return True


def nested_dict_get(pdict, keys, default=None):
    """
    Retrieve a value from a nested dict using a tuple of keys.

    If a is a dict, and a['b'] = {'c': 'd'},
    then nested_dict_get(a, ('b', 'c')) returns 'd'.
    """

    val = pdict
    for key in always_iterable(keys):
        try:
            val = val[key]
        except KeyError:
            return default
    return val


def get_listed_subparam(pdict, parent_param, subparam, default=None):
    """
    Returns nested_dict_get(pdict, (parent_param,subparam), default) if
    subparam is an entry in nested_dict_get(pdict, (parent_param, 'list'), [])

    This is a common idiom in Enzo-E's parameter parsing
    """
    if subparam in nested_dict_get(pdict, (parent_param, "list"), []):
        return nested_dict_get(pdict, (parent_param, subparam), default)
    return default


def get_particle_mass_correction(ds):
    """
    Normalize particle masses by the root grid cell volume.

    This correction is used for Enzo-E datasets where particle
    masses are stored as densities.
    """

    return (ds.domain_width / ds.domain_dimensions).prod() / ds.length_unit**3

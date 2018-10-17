import numpy as np
from walker import main, bundle_grids

def test_walker_donothing():
    level_seq = [    1,     1,     1,     1,     1,     1]
    ix_seq    = [(1,1), (2,1), (1,2), (2,2), (1,3), (2,3)]
    nblocks = np.array([2,3])
    patches = bundle_grids(levels=level_seq, ixs=ix_seq)
    res = main(patches, nblocks, ngrids=6)
    expected = patches
    assert deep_equality(res, expected)

def test_walker_missing_block():
    level_seq = [    2,     2,     2,     2]
    ix_seq    = [(1,1), (2,1), (1,2), (2,2)]
    nblocks = np.array([2,2])
    patches = bundle_grids(level_seq, ix_seq)
    res = main(patches, nblocks, ngrids=5)
    expected = bundle_grids(
        levels=[    1,     2,     2,     2,     2],
        ixs=   [(1,1), (1,1), (2,1), (1,2), (2,2)]
    )
    assert deep_equality(res, expected)

def test_walker_high_levels():
    level_seq = [    1,     2,     3,     4,     4,     4,     4,     2]
    ix_seq    = [(1,1), (2,1), (4,1), (7,3), (8,3), (7,4), (8,4), (2,2)]
    nblocks = np.array([1,1])
    patches = bundle_grids(level_seq, ix_seq)
    res = main(patches, nblocks, ngrids=9)
    expected = bundle_grids(
        levels=[    1,     2,     3,     3,     4,     4,     4,     4,     2],
        ixs=   [(1,1), (2,1), (4,1), (4,2), (7,3), (8,3), (7,4), (8,4), (2,2)]
    )
    for r,e in zip(res, expected):
        print(e)
        print(r)
        print()
    assert deep_equality(res, expected)

def deep_equality(result, expected):
    for r,e in zip(result, expected):
        for k in ('lvl', 'ix'):
            if not np.all(r[k] == e[k]):
                return False
    return True


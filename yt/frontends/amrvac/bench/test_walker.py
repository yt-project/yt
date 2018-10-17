import numpy as np
from walker import main, bundle_grids

def test_walker_simple():
    level_seq = [1] * 6
    ix_seq = [(1,1), (2,1), (1,2), (2,2), (1,3), (2,3)]
    nblocks = np.array([2,3])
    patches = bundle_grids(levels=level_seq, ixs=ix_seq)
    res = main(patches, nblocks, ngrids=6)
    expected = patches
    assert deep_equality(res, expected)

def test_walker_skip_levels():
    level_seq = [    2,     2,     2,     3,     3,     3,     3]
    ix_seq    = [(1,1), (1,2), (2,1), (3,3), (3,4), (4,3), (4,4)]
    nblocks = np.array([2,2])
    patches = bundle_grids(level_seq, ix_seq)
    res = main(patches, nblocks, ngrids=9)
    expected = bundle_grids(
        levels=[1,2,2,2,2,3,3,3,3],
        ixs=[(1,1), (1,1), (1,2), (2,1), (2,2), (3,3), (3,4), (4,3), (4,4)]
    )
    assert deep_equality(res, expected)

def deep_equality(result, expected):
    for r,e in zip(result, expected):
        for k in ('lvl', 'ix'):
            if not np.all(r[k] == e[k]):
                return False
    return True


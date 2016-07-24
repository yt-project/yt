from yt.testing import \
    fake_random_ds, \
    assert_equal
from yt.units.yt_array import \
    uconcatenate

def _get_dobjs(c):
    dobjs = [("sphere", ("center", (1.0, "unitary"))),
             ("sphere", ("center", (0.1, "unitary"))),
             ("ortho_ray", (0, (c[1], c[2]))),
             ("slice", (0, c[0])),
             #("disk", ("center", [0.1, 0.3, 0.6],
             #           (0.2, 'unitary'), (0.1, 'unitary'))),
             ("cutting", ([0.1, 0.3, 0.6], 'center')),
             ("all_data", ()),
            ]
    return dobjs

def test_chunking():
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs = nprocs)
        c = (ds.domain_right_edge + ds.domain_left_edge)/2.0 
        c += ds.arr(0.5/ds.domain_dimensions, "code_length")
        for dobj in _get_dobjs(c):
            obj = getattr(ds, dobj[0])(*dobj[1])
            coords = {'f':{}, 'i':{}}
            for t in ["io", "all", "spatial"]:
                coords['i'][t] = []
                coords['f'][t] = []
                for chunk in obj.chunks(None, t):
                    coords['f'][t].append(chunk.fcoords[:,:])
                    coords['i'][t].append(chunk.icoords[:,:])
                coords['f'][t] = uconcatenate(coords['f'][t])
                coords['i'][t] = uconcatenate(coords['i'][t])
                coords['f'][t].sort()
                coords['i'][t].sort()
            yield assert_equal, coords['f']['io'], coords['f']['all']
            yield assert_equal, coords['f']['io'], coords['f']['spatial']
            yield assert_equal, coords['i']['io'], coords['i']['all']
            yield assert_equal, coords['i']['io'], coords['i']['spatial']

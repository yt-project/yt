from yt.testing import *

def _get_dobjs(c):
    dobjs = [("sphere", ("center", (1.0, "unitary"))),
             ("sphere", ("center", (0.1, "unitary"))),
             ("ortho_ray", (0, (c[x_dict[0]], c[y_dict[0]]))),
             ("slice", (0, c[0])),
             #("disk", ("center", [0.1, 0.3, 0.6],
             #           (0.2, 'unitary'), (0.1, 'unitary'))),
             ("cutting", ([0.1, 0.3, 0.6], 'center')),
             ("all_data", ()),
            ]
    return dobjs

def test_chunking():
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs = nprocs)
        c = (pf.domain_right_edge + pf.domain_left_edge)/2.0 
        c += 0.5/pf.domain_dimensions
        for dobj in _get_dobjs(c):
            obj = getattr(pf.h, dobj[0])(*dobj[1])
            coords = {'f':{}, 'i':{}}
            for t in ["io", "all", "spatial"]:
                coords['i'][t] = []
                coords['f'][t] = []
                for chunk in obj.chunks(None, t):
                    coords['f'][t].append(chunk.fcoords[:,:])
                    coords['i'][t].append(chunk.icoords[:,:])
                coords['f'][t] = np.concatenate(coords['f'][t])
                coords['i'][t] = np.concatenate(coords['i'][t])
                coords['f'][t].sort()
                coords['i'][t].sort()
            yield assert_equal, coords['f']['io'], coords['f']['all']
            yield assert_equal, coords['f']['io'], coords['f']['spatial']
            yield assert_equal, coords['i']['io'], coords['i']['all']
            yield assert_equal, coords['i']['io'], coords['i']['spatial']

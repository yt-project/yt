import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_array_less

def setup():
    from yt.config import ytcfg
    ytcfg["yt","loglevel"] = "50"
    ytcfg["yt","__withintesting"] = "True"

def _difference(x1, x2, dw):
    rel = x1 - x2
    rel[rel >  dw/2.0] -= dw
    rel[rel < -dw/2.0] += dw
    return rel

def test_ellipsoid():
    # We decompose in different ways
    cs = [
          np.array([0.5, 0.5, 0.5]),
          np.array([0.1, 0.2, 0.3]),
          np.array([0.8, 0.8, 0.8])
          ]
    np.random.seed(int(0x4d3d3d3))
    for nprocs in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs = nprocs)
        DW = ds.domain_right_edge - ds.domain_left_edge
        min_dx = 2.0/ds.domain_dimensions
        ABC = np.random.random((3, 12)) * 0.1
        e0s = np.random.random((3, 12))
        tilts = np.random.random(12)
        ABC[:,0] = 0.1
        for i in range(12):
            for c in cs:
                A, B, C = reversed(sorted(ABC[:,i]))
                A = max(A, min_dx[0])
                B = max(B, min_dx[1])
                C = max(C, min_dx[2])
                e0 = e0s[:,i]
                tilt = tilts[i]
                ell = ds.ellipsoid(c, A, B, C, e0, tilt)
                yield assert_array_less, ell["radius"], A
                p = np.array([ell[ax] for ax in 'xyz'])
                dot_evec = [np.zeros_like(ell["radius"]) for i in range(3)]
                vecs = [ell._e0, ell._e1, ell._e2]
                mags = [ell._A, ell._B, ell._C]
                my_c = np.array([c]*p.shape[1]).transpose()
                dot_evec = [de.to_ndarray() for de in dot_evec]
                mags = [m.to_ndarray() for m in mags]
                for ax_i in range(3):
                    dist = _difference(p[ax_i,:], my_c[ax_i,:], DW[ax_i])
                    for ax_j in range(3):
                        dot_evec[ax_j] += dist * vecs[ax_j][ax_i]
                dist = 0
                for ax_i in range(3):
                    dist += dot_evec[ax_i]**2.0 / mags[ax_i]**2.0
                yield assert_array_less, dist, 1.0

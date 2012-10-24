from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

def test_ellipsoid():
    # We decompose in different ways
    cs = [[0.5, 0.5, 0.5],
          [0.1, 0.2, 0.3],
          [0.8, 0.8, 0.8]]
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs = nprocs)
        ABC = np.random.random((3, 12)) * 0.1
        e0s = np.random.random((3, 12))
        tilts = np.random.random((3, 12))
        ABC[:,0] = 0.1

        for i in range(12):
            for c in cs:
                A, B, C = reversed(sorted(ABC[:,i]))
                e0 = e0s[:,i]
                tilt = tilts[:,i]
                ell = pf.h.ellipsoid(c, A, B, C, e0, tilt)
                yield assert_equal, np.all(ell["Radius"] <= A), True
                for i in xrange(ell["Radius"].size):
                    pos = np.array([ell[ax][i] for ax in 'xyz'])
                    v = 0.0
                    v += (pos * ell._e0).sum() / ell._A
                    v += (pos * ell._e1).sum() / ell._B
                    v += (pos * ell._e2).sum() / ell._C
                    yield assert_equal, (v <= 1.0), True

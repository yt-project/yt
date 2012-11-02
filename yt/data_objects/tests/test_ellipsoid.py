from yt.testing import *

def setup():
    from yt.config import ytcfg
    ytcfg["yt","loglevel"] = "50"
    ytcfg["yt","__withintesting"] = "True"

def test_ellipsoid():
    # We decompose in different ways
    cs = [np.array([0.5, 0.5, 0.5]),
          np.array([0.1, 0.2, 0.3]),
          np.array([0.8, 0.8, 0.8])]
    for nprocs in [1, 2, 4, 8]:
        pf = fake_random_pf(64, nprocs = nprocs)
        min_dx = 2.0/pf.domain_dimensions
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
                ell = pf.h.ellipsoid(c, A, B, C, e0, tilt)
                yield assert_equal, np.all(ell["Radius"] <= A), True
                p = np.array([ell[ax] for ax in 'xyz'])
                v  = np.zeros_like(ell["Radius"])
                v += (((p - c[:,None]) * ell._e0[:,None]).sum(axis=0) / ell._A)**2
                v += (((p - c[:,None]) * ell._e1[:,None]).sum(axis=0) / ell._B)**2
                v += (((p - c[:,None]) * ell._e2[:,None]).sum(axis=0) / ell._C)**2
                yield assert_equal, np.all(np.sqrt(v) <= 1.0), True

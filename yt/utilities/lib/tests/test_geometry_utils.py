import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_array_less, \
    assert_array_equal
from yt.utilities.lib.misc_utilities import obtain_rvec, obtain_rv_vec

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")

def test_get_morton_indices():
    from yt.utilities.lib.geometry_utils import get_morton_indices,get_morton_indices_unravel
    INDEX_MAX_64 = np.uint64(2097151)
    li = np.arange(6,dtype=np.uint64).reshape((2,3))
    mi_ans = np.array([10,229],dtype=np.uint64)
    mi_out = get_morton_indices(li)
    mi_out2 = get_morton_indices_unravel(li[:,0],li[:,1],li[:,2])
    assert_array_equal(mi_out,mi_ans)
    assert_array_equal(mi_out2,mi_ans)
    li[0,:] = INDEX_MAX_64*np.ones(3,dtype=np.uint64)
    assert_raises(ValueError,get_morton_indices,li)
    assert_raises(ValueError,get_morton_indices_unravel,li[:,0],li[:,1],li[:,2])

def test_get_morton_points():
    from yt.utilities.lib.geometry_utils import get_morton_points
    mi = np.array([10,229],dtype=np.uint64)
    li_ans = np.arange(6,dtype=np.uint64).reshape((2,3))
    li_out = get_morton_points(mi)
    assert_array_equal(li_out,li_ans)

def test_compare_morton():
    # TODO: Add error messages to assertions
    from yt.utilities.lib.geometry_utils import compare_morton
    # Diagonal
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([1.0,1.0,1.0],dtype=np.float64)
    assert_equal(compare_morton(p,q),1)
    assert_equal(compare_morton(q,p),0)
    assert_equal(compare_morton(p,p),0)
    # 1-1 vs 0-1
    p = np.array([1.0,1.0,0.0],dtype=np.float64)
    q = np.array([1.0,1.0,1.0],dtype=np.float64)
    assert_equal(compare_morton(p,q),1)
    assert_equal(compare_morton(q,p),0)
    assert_equal(compare_morton(p,p),0)
    # x advance, y decrease
    p = np.array([0.0,1.0,0.0],dtype=np.float64)
    q = np.array([1.0,0.0,0.0],dtype=np.float64)
    assert_equal(compare_morton(p,q),1)
    assert_equal(compare_morton(q,p),0)
    assert_equal(compare_morton(p,p),0)
    # x&y advance, z decrease
    p = np.array([0.0,0.0,1.0],dtype=np.float64)
    q = np.array([1.0,1.0,0.0],dtype=np.float64)
    assert_equal(compare_morton(p,q),1)
    assert_equal(compare_morton(q,p),0)
    assert_equal(compare_morton(p,p),0)

def test_get_morton_argsort():
    from yt.utilities.lib.geometry_utils import get_morton_argsort,get_morton_points
    N = 10
    mi = np.hstack((np.arange(N,dtype=np.uint64),np.array([646904,780276],dtype=np.uint64)))
    N = len(mi)
    pos = get_morton_points(mi).astype(np.float64)
    sort_out = np.arange(N,dtype=np.uint64)
    sort_shf = np.arange(N,dtype=np.uint64)
    np.random.shuffle(sort_shf)
    sort_ans = np.argsort(sort_shf)
    get_morton_argsort(pos[sort_shf,:],0,N-1,sort_out)
    assert_array_equal(sort_out,sort_ans)

def test_dist():
    from yt.utilities.lib.geometry_utils import dist
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([0.0,0.0,0.0],dtype=np.float64)
    assert_equal(dist(p,q),0.0)
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([1.0,0.0,0.0],dtype=np.float64)
    assert_equal(dist(p,q),1.0)
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([1.0,1.0,0.0],dtype=np.float64)
    assert_equal(dist(p,q),np.sqrt(2.0))
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([1.0,1.0,1.0],dtype=np.float64)
    assert_equal(dist(p,q),np.sqrt(3.0))

def test_knn_direct():
    from yt.utilities.lib.geometry_utils import knn_direct
    k = 5
    N = 2*k
    idx = np.arange(N,dtype=np.uint64)
    rad = np.arange(N,dtype=np.float64)
    pos = np.vstack(3*[rad**2/3.0]).T
    sort_shf = np.arange(N,dtype=np.uint64)
    np.random.shuffle(sort_shf)
    sort_ans = np.argsort(sort_shf)[:k]
    sort_out = knn_direct(pos[sort_shf,:], k, sort_ans[0], idx)
    assert_array_equal(sort_out,sort_ans)

# TODO: test of quadtree (.pxd)

def point_grid(n_per_dim,ndim):
    q = np.arange(n_per_dim,dtype=np.float64)+0.5 # Middle of each cell
    out = np.meshgrid(*tuple(ndim*[q]))
    pos = np.vstack(tuple([iout.flatten() for iout in out])).T
    DLE = np.array(ndim*[0.0],dtype=np.float64)
    DRE = np.array(ndim*[n_per_dim],dtype=np.float64)
    pos = pos.astype(np.float64)
    # ind = np.arange(pos.shape[0])
    # np.random.shuffle(ind)
    # pos = pos[ind,:]
    return pos,DLE,DRE

def point_random(n_per_dim,ndim,seed=1):
    np.random.seed(seed)
    pos = np.random.random_sample((n_per_dim**ndim,ndim)).astype(np.float64)
    DLE = np.array(ndim*[0.0],dtype=np.float64)
    DRE = np.array(ndim*[1.0],dtype=np.float64)
    pos = pos.astype(np.float64)
    return pos,DLE,DRE

def test_knn_morton():
    from yt.utilities.lib.geometry_utils import knn_direct,knn_morton,get_morton_argsort,dist
    Np_d = 50 # 100
    Nd = 3
    Np = Np_d**Nd
    idx_test = np.uint64(Np/2)
    k = 6 # 27
    idx_notest = np.hstack([np.arange(idx_test),np.arange((idx_test+1),Np)]).astype(np.uint64)
    # Random points
    pos,DLE,DRE = point_random(Np_d,Nd)
    sort_fwd = np.arange(pos.shape[0],dtype=np.uint64)
    get_morton_argsort(pos,0,pos.shape[0]-1,sort_fwd)
    sort_rev = np.argsort(sort_fwd).astype(np.uint64)
    knn_dir = sorted(knn_direct(pos[sort_fwd,:],k,sort_rev[idx_test],sort_rev[idx_notest]))
    knn_mor = sorted(knn_morton(pos[sort_fwd,:],k,sort_rev[idx_test],DLE=DLE,DRE=DRE,issorted=True))
    if True:
        print pos[idx_test,:]
        id1 = 14252
        id2 = 14502
        print id1,pos[sort_fwd[id1],:],dist(pos[idx_test,:],pos[sort_fwd[id1],:])
        print id2,pos[sort_fwd[id2],:],dist(pos[idx_test,:],pos[sort_fwd[id2],:])
    assert_array_equal(knn_mor,knn_dir,
                       err_msg="{} random points & k = {}".format(Np,k))
    # Grid points
    pos,DLE,DRE = point_grid(Np_d,Nd)
    sort_fwd = np.arange(pos.shape[0],dtype=np.uint64)
    get_morton_argsort(pos,0,pos.shape[0]-1,sort_fwd)
    sort_rev = np.argsort(sort_fwd).astype(np.uint64)
    knn_dir = sorted(knn_direct(pos[sort_fwd,:],k,sort_rev[idx_test],sort_rev[idx_notest]))
    knn_mor = sorted(knn_morton(pos[sort_fwd,:],k,sort_rev[idx_test],DLE=DLE,DRE=DRE,issorted=True))
    if True:
        print pos[idx_test,:]
        id1 = 14252
        id2 = 14502
        print id1,pos[sort_fwd[id1],:],dist(pos[idx_test,:],pos[sort_fwd[id1],:])
        print id2,pos[sort_fwd[id2],:],dist(pos[idx_test,:],pos[sort_fwd[id2],:])
    assert_array_equal(knn_mor,knn_dir,
                       err_msg="grid of {} points & k = {}".format(Np,k))


def test_csearch_morton():
    from yt.utilities.lib.geometry_utils import get_morton_argsort,csearch_morton,ORDER_MAX
    k = 5
    i = 0
    xN = 3
    N = xN**3
    xf = np.arange(xN,dtype=np.float64)+1
    DLE = np.zeros(3,dtype=np.float64)
    DRE = (xN+1)*np.ones(3,dtype=np.float64)
    X,Y,Z = np.meshgrid(xf,xf,xf)
    pos = np.vstack((X.flatten(),Y.flatten(),Z.flatten())).T
    sort_fwd = np.arange(N, dtype=np.uint64)
    get_morton_argsort(pos,0,N-1,sort_fwd)
    out = csearch_morton(pos[sort_fwd,:], k, i, Ai, l, h, order, DLE, DRE, ORDER_MAX)
    # assert_array_equal(out,ans)

def test_obtain_rvec():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])
    
    dd = ds.sphere((0.5,0.5,0.5), 0.2)

    coords = obtain_rvec(dd)

    r = np.sqrt(np.sum(coords*coords,axis=0))

    assert_array_less(r.max(), 0.2)

    assert_array_less(0.0, r.min())

def test_obtain_rv_vec():
    ds = fake_random_ds(64, nprocs=8, fields=_fields, 
           negative = [False, True, True, True])

    dd = ds.all_data()

    vels = obtain_rv_vec(dd)

    assert_array_equal(vels[0,:], dd['velocity_x'])
    assert_array_equal(vels[1,:], dd['velocity_y'])
    assert_array_equal(vels[2,:], dd['velocity_z'])

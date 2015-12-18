from yt.testing import *
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
    from yt.utilities.lib.geometry_utils import compare_morton
    p = np.array([0.0,0.0,0.0],dtype=np.float64)
    q = np.array([1.0,1.0,1.0],dtype=np.float64)
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
    print pos[sort_shf,:]
    get_morton_argsort(pos[sort_shf,:],0,N-1,sort_out)
    assert_array_equal(sort_out,sort_ans)

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

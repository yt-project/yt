from yt.testing import *
from yt.utilities.lib.misc_utilities import obtain_rvec, obtain_rv_vec

_fields = ("density", "velocity_x", "velocity_y", "velocity_z")

# TODO: error compact/spread bits for incorrect size
# TODO: test msdb for [0,0], [1,1], [2,2] etc.

def test_spread_bits():
    from yt.utilities.lib.geometry_utils import spread_bits
    li = [(np.uint64(0b111111111111111111111), np.uint64(0b1001001001001001001001001001001001001001001001001001001001001))]
    for i,ans in li:
        out = spread_bits(i)
        assert_equal(out,ans)

def test_compact_bits():
    from yt.utilities.lib.geometry_utils import compact_bits
    li = [(np.uint64(0b111111111111111111111), np.uint64(0b1001001001001001001001001001001001001001001001001001001001001))]
    for ans,i in li:
        out = compact_bits(i)
        assert_equal(out,ans)

def test_spread_and_compact_bits():
    from yt.utilities.lib.geometry_utils import spread_bits,compact_bits
    li = [np.uint64(0b111111111111111111111)]
    for ans in li:
        mi = spread_bits(ans)
        out = compact_bits(mi)
        assert_equal(out,ans)

def test_lsz():
    from yt.utilities.lib.geometry_utils import lsz
    li = [(np.uint64(0b1001001001001001001001001001001001001001001001001001001001001)  ,3*21, 3, 0),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001001001000)  , 3*0, 3, 0),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001001000001)  , 3*1, 3, 0),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001000001001)  , 3*2, 3, 0),
          (np.uint64(0b10010010010010010010010010010010010010010010010010010010010010) , 3*0, 3, 0),
          (np.uint64(0b100100100100100100100100100100100100100100100100100100100100100), 3*0, 3, 0),
          (np.uint64(0b100), 0, 1, 0),
          (np.uint64(0b100), 1, 1, 1),
          (np.uint64(0b100), 3, 1, 2),
          (np.uint64(0b100), 3, 1, 3)]
    for i,ans,stride,start in li:
        out = lsz(i,stride=stride,start=start)
        assert_equal(out,ans)

def test_lsb():
    from yt.utilities.lib.geometry_utils import lsb
    li = [(np.uint64(0b1001001001001001001001001001001001001001001001001001001001001)  , 3*0),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001001001000)  , 3*1),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001001000000)  , 3*2),
          (np.uint64(0b1001001001001001001001001001001001001001001001001001000000000)  , 3*3),
          (np.uint64(0b10010010010010010010010010010010010010010010010010010010010010) ,3*21),
          (np.uint64(0b100100100100100100100100100100100100100100100100100100100100100),3*21)]
    for i,ans in li:
        out = lsb(i,stride=3)
        assert_equal(out,ans)

def test_bitwise_addition():
    from yt.utilities.lib.geometry_utils import bitwise_addition
    # TODO: Handle negative & periodic boundaries
    begin = 1
    end = 5
    lz = [(0,1),
#          (0,-1),
          (1,1),
          (1,2),
          (1,4),
          (1,-1),
          (2,1),
          (2,2),
          (2,-1),
          (2,-2),
          (3,1),
          (3,5),
          (3,-1)]
    for i,a in lz:
        i = np.uint64(i)
        a = np.int64(a)
        out = bitwise_addition(i,a,stride=1,start=0)
        # print bin(i),bin(a),bin(np.uint64(i+a)),bin(out)
        assert_equal(out,i+a)

#def test_add_to_morton_coord():
#    from yt.utilities.lib.geometry_utils import add_to_morton_coord
    

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

def test_morton_qsort(seed=1,recursive=False,use_loop=False):
    from yt.utilities.lib.geometry_utils import morton_qsort,get_morton_points
    np.random.seed(seed)
    N = 100
    mi = np.hstack((np.arange(N,dtype=np.uint64),#+600000,
                    np.array([646904,780276],dtype=np.uint64)))
    N = len(mi)
    pos = get_morton_points(mi).astype(np.float64)
    sort_out = np.arange(N,dtype=np.uint64)
    sort_shf = np.arange(N,dtype=np.uint64)
    np.random.shuffle(sort_shf)
    sort_ans = np.argsort(sort_shf)
    morton_qsort(pos[sort_shf,:],0,N-1,sort_out,
                 recursive=recursive,use_loop=use_loop)
    if True:
        from yt.utilities.lib.geometry_utils import compare_morton,xor_msb_cy,msdb_cy,ifrexp_cy
        count = 0
        for ind1 in range(N):
            if sort_out[ind1]!=sort_ans[ind1]:
                for ind2 in range(ind1+1,N):
                    i1 = sort_ans[ind1]
                    i2 = sort_ans[ind2]
                    p1 = pos[sort_shf[i1],:]
                    p2 = pos[sort_shf[i2],:]
                    if not compare_morton(p1,p2):
                        count+=1
                        print compare_morton(p1,p2),compare_morton(p2,p1)
                        print 'Found {}'.format(count)
                        print '    ',i1, mi[sort_shf[i1]], p1
                        print '    ',i2, mi[sort_shf[i2]], p2
                        for j in range(3):
                            im1,ie1 = ifrexp_cy(p1[j])
                            im2,ie2 = ifrexp_cy(p2[j])
                            print '    ',j,xor_msb_cy(p1[j],p2[j]),[ie1,ie2],[im1,im2],msdb_cy(im1,im2)
    assert_array_equal(sort_out,sort_ans)

def test_morton_neighbor():
    from yt.utilities.lib.geometry_utils import morton_neighbor, get_morton_indices, get_morton_index
    order = 20
    imax = 1 << order
    p = np.array([[imax/2,imax/2,imax/2],
                  [imax/2,imax/2,0     ],
                  [imax/2,imax/2,imax  ]],dtype=np.uint64)
    add = np.array([[+1, 0, 0],
                    [+1,+1, 0],[+1,+1,+1],[+1,+1,-1],
                    [+1,-1, 0],[+1,-1,+1],[+1,-1,-1],
                    [+1, 0,+1],[+1, 0,-1],
                    [-1, 0, 0],
                    [-1,+1, 0],[-1,+1,+1],[-1,+1,-1],
                    [-1,-1, 0],[-1,-1,+1],[-1,-1,-1],
                    [-1, 0,+1],[-1, 0,-1],
                    [ 0,+1, 0],
                    [ 0,+1,+1],[ 0,+1,-1],
                    [ 0,-1, 0],
                    [ 0,-1,+1],[ 0,-1,-1],
                    [ 0, 0,+1],
                    [ 0, 0,-1]])
    p_ans = np.array([[imax/2,imax/2,imax/2+1],
                      [imax/2,imax/2,imax/2-1],
                      [imax/2,imax/2,imax-1  ],
                      [imax/2,imax/2,1       ],
                      [imax/2,imax/2+1,imax/2+1],
                      [imax/2-1,imax/2-1,imax/2],
                      [imax/2-1,imax/2,imax/2+1],
                      [imax/2,imax/2-1,imax-1  ],
                      [imax/2,imax/2+1,1       ]],dtype=np.uint64)
    mi_ans = get_morton_indices(p_ans)
    assert_equal(morton_neighbor(p[0,:],[2],[+1],order=order),mi_ans[0])
    assert_equal(morton_neighbor(p[0,:],[2],[-1],order=order),mi_ans[1])
    assert_equal(morton_neighbor(p[1,:],[2],[-1],order=order,periodic=False),-1)
    assert_equal(morton_neighbor(p[2,:],[2],[+1],order=order,periodic=False),-1)
    assert_equal(morton_neighbor(p[1,:],[2],[-1],order=order,periodic=True ),mi_ans[2])
    assert_equal(morton_neighbor(p[2,:],[2],[+1],order=order,periodic=True ),mi_ans[3])
    assert_equal(morton_neighbor(p[0,:],[1,2],[+1,+1],order=order),mi_ans[4])
    assert_equal(morton_neighbor(p[0,:],[0,1],[-1,-1],order=order),mi_ans[5])
    assert_equal(morton_neighbor(p[0,:],[0,2],[-1,+1],order=order),mi_ans[6])
    assert_equal(morton_neighbor(p[1,:],[1,2],[-1,-1],order=order,periodic=False),-1)
    assert_equal(morton_neighbor(p[2,:],[1,2],[+1,+1],order=order,periodic=False),-1)
    assert_equal(morton_neighbor(p[1,:],[1,2],[-1,-1],order=order,periodic=True ),mi_ans[7])
    assert_equal(morton_neighbor(p[2,:],[1,2],[+1,+1],order=order,periodic=True ),mi_ans[8])

def test_get_morton_neighbors():
    from yt.utilities.lib.geometry_utils import get_morton_neighbors, get_morton_indices
    order = 20
    imax = 1 << order
    p = np.array([[imax/2,imax/2,imax/2],
                  [imax/2,imax/2,0     ],
                  [imax/2,imax/2,imax  ]],dtype=np.uint64)
    pn_non = [
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,imax/2],
            [imax/2+1,imax/2+1,imax/2],[imax/2+1,imax/2+1,imax/2+1],[imax/2+1,imax/2+1,imax/2-1],
            [imax/2+1,imax/2-1,imax/2],[imax/2+1,imax/2-1,imax/2+1],[imax/2+1,imax/2-1,imax/2-1],
            [imax/2+1,imax/2,imax/2+1],[imax/2+1,imax/2,imax/2-1],
            [imax/2-1,imax/2,imax/2],
            [imax/2-1,imax/2+1,imax/2],[imax/2-1,imax/2+1,imax/2+1],[imax/2-1,imax/2+1,imax/2-1],
            [imax/2-1,imax/2-1,imax/2],[imax/2-1,imax/2-1,imax/2+1],[imax/2-1,imax/2-1,imax/2-1],
            [imax/2-1,imax/2,imax/2+1],[imax/2-1,imax/2,imax/2-1],
            # y +/- 1
            [imax/2,imax/2+1,imax/2],
            [imax/2,imax/2+1,imax/2+1],[imax/2,imax/2+1,imax/2-1],
            [imax/2,imax/2-1,imax/2],
            [imax/2,imax/2-1,imax/2+1],[imax/2,imax/2-1,imax/2-1],
            # x +/- 1
            [imax/2,imax/2,imax/2+1],
            [imax/2,imax/2,imax/2-1]],dtype=np.uint64),
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,0],
            [imax/2+1,imax/2+1,0],[imax/2+1,imax/2+1,1],
            [imax/2+1,imax/2-1,0],[imax/2+1,imax/2-1,1],
            [imax/2+1,imax/2,1],
            [imax/2-1,imax/2,0],
            [imax/2-1,imax/2+1,0],[imax/2-1,imax/2+1,1],
            [imax/2-1,imax/2-1,0],[imax/2-1,imax/2-1,1],
            [imax/2-1,imax/2,1],
            # y +/- 1
            [imax/2,imax/2+1,0],
            [imax/2,imax/2+1,1],
            [imax/2,imax/2-1,0],
            [imax/2,imax/2-1,1],
            # z +/- 1
            [imax/2,imax/2,0+1]],dtype=np.uint64),
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,imax],
            [imax/2+1,imax/2+1,imax],[imax/2+1,imax/2+1,imax-1],
            [imax/2+1,imax/2-1,imax],[imax/2+1,imax/2-1,imax-1],
            [imax/2+1,imax/2,imax-1],
            [imax/2-1,imax/2,imax],
            [imax/2-1,imax/2+1,imax],[imax/2-1,imax/2+1,imax-1],
            [imax/2-1,imax/2-1,imax],[imax/2-1,imax/2-1,imax-1],
            [imax/2-1,imax/2,imax-1],
            # y +/- 1
            [imax/2,imax/2+1,imax],
            [imax/2,imax/2+1,imax-1],
            [imax/2,imax/2-1,imax],
            [imax/2,imax/2-1,imax-1],
            # z +/- 1
            [imax/2,imax/2,imax-1]],dtype=np.uint64)]
    pn_per = [
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,imax/2],
            [imax/2+1,imax/2+1,imax/2],[imax/2+1,imax/2+1,imax/2+1],[imax/2+1,imax/2+1,imax/2-1],
            [imax/2+1,imax/2-1,imax/2],[imax/2+1,imax/2-1,imax/2+1],[imax/2+1,imax/2-1,imax/2-1],
            [imax/2+1,imax/2,imax/2+1],[imax/2+1,imax/2,imax/2-1],
            [imax/2-1,imax/2,imax/2],
            [imax/2-1,imax/2+1,imax/2],[imax/2-1,imax/2+1,imax/2+1],[imax/2-1,imax/2+1,imax/2-1],
            [imax/2-1,imax/2-1,imax/2],[imax/2-1,imax/2-1,imax/2+1],[imax/2-1,imax/2-1,imax/2-1],
            [imax/2-1,imax/2,imax/2+1],[imax/2-1,imax/2,imax/2-1],
            # y +/- 1
            [imax/2,imax/2+1,imax/2],
            [imax/2,imax/2+1,imax/2+1],[imax/2,imax/2+1,imax/2-1],
            [imax/2,imax/2-1,imax/2],
            [imax/2,imax/2-1,imax/2+1],[imax/2,imax/2-1,imax/2-1],
            # z +/- 1
            [imax/2,imax/2,imax/2+1],
            [imax/2,imax/2,imax/2-1]],dtype=np.uint64),
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,0],
            [imax/2+1,imax/2+1,0],[imax/2+1,imax/2+1,1],[imax/2+1,imax/2+1,imax-1],
            [imax/2+1,imax/2-1,0],[imax/2+1,imax/2-1,1],[imax/2+1,imax/2-1,imax-1],
            [imax/2+1,imax/2,1],[imax/2+1,imax/2,imax-1],
            [imax/2-1,imax/2,0],
            [imax/2-1,imax/2+1,0],[imax/2-1,imax/2+1,1],[imax/2-1,imax/2+1,imax-1],
            [imax/2-1,imax/2-1,0],[imax/2-1,imax/2-1,1],[imax/2-1,imax/2-1,imax-1],
            [imax/2-1,imax/2,1],[imax/2-1,imax/2,imax-1],
            # y +/- 1
            [imax/2,imax/2+1,0],
            [imax/2,imax/2+1,1],[imax/2,imax/2+1,imax-1],
            [imax/2,imax/2-1,0],
            [imax/2,imax/2-1,1],[imax/2,imax/2-1,imax-1],
            # z +/- 1
            [imax/2,imax/2,0+1],
            [imax/2,imax/2,imax-1]],dtype=np.uint64),
        np.array([
            # x +/- 1
            [imax/2+1,imax/2,imax],
            [imax/2+1,imax/2+1,imax],[imax/2+1,imax/2+1,1],[imax/2+1,imax/2+1,imax-1],
            [imax/2+1,imax/2-1,imax],[imax/2+1,imax/2-1,1],[imax/2+1,imax/2-1,imax-1],
            [imax/2+1,imax/2,1],[imax/2+1,imax/2,imax-1],
            [imax/2-1,imax/2,imax],
            [imax/2-1,imax/2+1,imax],[imax/2-1,imax/2+1,1],[imax/2-1,imax/2+1,imax-1],
            [imax/2-1,imax/2-1,imax],[imax/2-1,imax/2-1,1],[imax/2-1,imax/2-1,imax-1],
            [imax/2-1,imax/2,1],[imax/2-1,imax/2,imax-1],
            # y +/- 1
            [imax/2,imax/2+1,imax],
            [imax/2,imax/2+1,1],[imax/2,imax/2+1,imax-1],
            [imax/2,imax/2-1,imax],
            [imax/2,imax/2-1,1],[imax/2,imax/2-1,imax-1],
            # z +/- 1
            [imax/2,imax/2,1],
            [imax/2,imax/2,imax-1]],dtype=np.uint64)]
    mi = get_morton_indices(p)
    N = mi.shape[0]
    # Non-periodic
    for i in range(N):
        out = get_morton_neighbors(np.array([mi[i]],dtype=np.uint64),order=order,periodic=False)
        ans = get_morton_indices(pn_non[i])
        assert_array_equal(out,ans,err_msg="Non-periodic: {}".format(i))
    # Periodic
    for i in range(N):
        out = get_morton_neighbors(np.array([mi[i]],dtype=np.uint64),order=order,periodic=True)
        ans = get_morton_indices(pn_per[i])
        assert_array_equal(out,ans,err_msg="Periodic: {}".format(i))

def time_bitwise_addition():
    import time
    from yt.utilities.lib.geometry_utils import get_morton_points, get_morton_indices, bitwise_addition
    xarr = np.array(range(100)+[np.uint64(0b11111111111111111111)],dtype=np.uint64)
    # Explicit spreading and compacting
    t1 = time.time()
    for x in xarr:
        p = get_morton_points(np.array([x],dtype=np.uint64))
        p[0]+=1
        p[0]+=2
        x1 = get_morton_indices(p)
    t2 = time.time()
    print("Explicit bit spreading/compacting: {:f}".format(t2-t1))
    # Bitwise addition
    t1 = time.time()
    for x in xarr:
        x2 = bitwise_addition(x,np.int64(1),stride=3,start=2)
        x2 = bitwise_addition(x,np.int64(2),stride=3,start=2)
    t2 = time.time()
    print("Using bitwise addition: {:f}".format(t2-t1))

def time_morton_qsort(seed=1):
    # Not the most effecient test, but not meant to be run much
    import time
    # Recursive, no loop
    t1 = time.time()
    test_morton_qsort(recursive=True,use_loop=False)
    t2 = time.time()
    print("Recursive, no loop: {:f}".format(t2-t1))
    # Recursive, loop
    t1 = time.time()
    test_morton_qsort(recursive=True,use_loop=True)
    t2 = time.time()
    print("Recursive, loop:    {:f}".format(t2-t1))
    # Iterative, no loop
    t1 = time.time()
    test_morton_qsort(recursive=False,use_loop=False)
    t2 = time.time()
    print("Iterative, no loop: {:f}".format(t2-t1))
    # Iterative, loop
    t1 = time.time()
    test_morton_qsort(recursive=False,use_loop=True)
    t2 = time.time()
    print("Iterative, loop:    {:f}".format(t2-t1))

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

def test_knn_direct(seed=1):
    from yt.utilities.lib.geometry_utils import knn_direct
    np.random.seed(seed)
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

def gen_points_grid(n_per_dim,ndim):
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

def gen_points_random(n_per_dim,ndim,seed=1):
    np.random.seed(seed)
    pos = np.random.random_sample((n_per_dim**ndim,ndim)).astype(np.float64)
    DLE = np.array(ndim*[0.0],dtype=np.float64)
    DRE = np.array(ndim*[1.0],dtype=np.float64)
    pos = pos.astype(np.float64)
    return pos,DLE,DRE

def gen_points(n_per_dim,ndim,seed=1,grid=False):
    if grid:
        return gen_points_grid(n_per_dim,ndim)
    else:
        return gen_points_random(n_per_dim,ndim,seed=seed)

def test_knn_morton():
    from yt.utilities.lib.geometry_utils import knn_direct,knn_morton,morton_qsort,dist
    Np_d = 10 # 100
    Nd = 3
    Np = Np_d**Nd
    idx_test = np.uint64(Np/2)
    k = 6 # 27
    idx_notest = np.hstack([np.arange(idx_test),np.arange((idx_test+1),Np)]).astype(np.uint64)
    # Random points
    pos,DLE,DRE = gen_points_random(Np_d,Nd)
    sort_fwd = np.arange(pos.shape[0],dtype=np.uint64)
    morton_qsort(pos,0,pos.shape[0]-1,sort_fwd)
    sort_rev = np.argsort(sort_fwd).astype(np.uint64)
    knn_dir = sorted(knn_direct(pos[sort_fwd,:],k,sort_rev[idx_test],sort_rev[idx_notest]))
    knn_mor = sorted(knn_morton(pos[sort_fwd,:],k,sort_rev[idx_test],DLE=DLE,DRE=DRE,issorted=True))
    if False:
        print pos[idx_test,:]
        id1 = 14252
        id2 = 14502
        print id1,pos[sort_fwd[id1],:],dist(pos[idx_test,:],pos[sort_fwd[id1],:])
        print id2,pos[sort_fwd[id2],:],dist(pos[idx_test,:],pos[sort_fwd[id2],:])
    assert_array_equal(knn_mor,knn_dir,
                       err_msg="{} random points & k = {}".format(Np,k))
    # Grid points
    pos,DLE,DRE = gen_points_grid(Np_d,Nd)
    sort_fwd = np.arange(pos.shape[0],dtype=np.uint64)
    morton_qsort(pos,0,pos.shape[0]-1,sort_fwd)
    sort_rev = np.argsort(sort_fwd).astype(np.uint64)
    knn_dir = sorted(knn_direct(pos[sort_fwd,:],k,sort_rev[idx_test],sort_rev[idx_notest]))
    knn_mor = sorted(knn_morton(pos[sort_fwd,:],k,sort_rev[idx_test],DLE=DLE,DRE=DRE,issorted=True))
    if False:
        print pos[idx_test,:]
        id1 = 14252
        id2 = 14502
        print id1,pos[sort_fwd[id1],:],dist(pos[idx_test,:],pos[sort_fwd[id1],:])
        print id2,pos[sort_fwd[id2],:],dist(pos[idx_test,:],pos[sort_fwd[id2],:])
    assert_array_equal(knn_mor,knn_dir,
                       err_msg="grid of {} points & k = {}".format(Np,k))


def test_csearch_morton():
    from yt.utilities.lib.geometry_utils import morton_qsort,csearch_morton,ORDER_MAX
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
    morton_qsort(pos,0,N-1,sort_fwd)
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

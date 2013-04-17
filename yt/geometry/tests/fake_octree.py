from yt.geometry.fake_octree import create_fake_octree
import numpy as np

max_leaf = 100
max_level = 5
dn = 2
dd = np.ones(3,dtype='i4')*dn
dle = np.ones(3,dtype='f8')*0.0
dre = np.ones(3,dtype='f8')
fsub = 0.90

octtree = create_fake_octree(max_leaf, max_level, dd, dle, dre, fsub)

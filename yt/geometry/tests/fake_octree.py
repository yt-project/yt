from yt.geometry.fake_octree import create_fake_octree
from yt.geometry.oct_container import RAMSESOctreeContainer
import numpy as np

nocts= 100
max_level = 12
dn = 2
dd = np.ones(3,dtype='i4')*dn
dle = np.ones(3,dtype='f8')*0.0
dre = np.ones(3,dtype='f8')
fsub = 0.25
domain = 1

oct_handler = RAMSESOctreeContainer(dd,dle,dre)
create_fake_octree(oct_handler, nocts, max_level, dd, dle, dre, fsub)
print "filled"
print oct_handler.check(1, print_all=1)
mask = np.ones((nocts,8),dtype='bool')
cell_count = nocts*8
level_counts = oct_handler.count_levels(max_level, 1, mask)
print level_counts
print "fcoords"
fc = oct_handler.fcoords(domain,mask,cell_count,level_counts)
print level_counts, level_counts.sum()
print [np.unique(fc[:,ax]).shape[0] for ax in range(3)]
print fc
print fc.shape
import pdb; pdb.set_trace()

#Now take the particles and recreate the same octree

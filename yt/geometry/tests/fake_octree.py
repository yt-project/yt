from yt.geometry.fake_octree import create_fake_octree
from yt.geometry.oct_container import RAMSESOctreeContainer
import numpy as np

nocts= 100
max_level = 12
dn = 2
dd = np.ones(3,dtype='i4')*dn
dle = np.ones(3,dtype='f8')*0.0
dre = np.ones(3,dtype='f8')
fsub = 0.10
domain = 0

oct_handler = RAMSESOctreeContainer(dd,dle,dre)
create_fake_octree(oct_handler, nocts, max_level, dd, dle, dre, fsub)
print "filled"
print oct_handler.check(domain, print_all=1)
mask = np.ones(nocts,dtype='bool')
print mask
cell_count = nocts*8
level_counts = np.array([nocts]) # not used anyway
fc = oct_handler.fcoords(domain,mask,cell_count)
print fc
print fc.shape


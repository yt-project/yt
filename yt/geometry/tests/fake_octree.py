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
mask = np.ones((nocts,8),dtype='bool')
cell_count = nocts*8
oct_counts = oct_handler.count_levels(max_level, 1, mask)
level_counts = np.concatenate(([0,],np.cumsum(oct_counts)))
fc = oct_handler.fcoords(domain,mask,cell_count, level_counts.copy())

#Now take the particles and recreate the same octree

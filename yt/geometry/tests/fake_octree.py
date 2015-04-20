from __future__ import print_function
from yt.geometry.fake_octree import create_fake_octree
from yt.geometry.oct_container import RAMSESOctreeContainer, ParticleOctreeContainer
import numpy as np

nocts= 3
max_level = 12
dn = 2
dd = np.ones(3,dtype='i4')*dn
dle = np.ones(3,dtype='f8')*0.0
dre = np.ones(3,dtype='f8')
fsub = 0.25
domain = 1

oct_handler = RAMSESOctreeContainer(dd,dle,dre)
leaves = create_fake_octree(oct_handler, nocts, max_level, dd, dle, dre, fsub)
mask = np.ones((nocts,8),dtype='bool')
cell_count = nocts*8
oct_counts = oct_handler.count_levels(max_level, 1, mask)
level_counts = np.concatenate(([0,],np.cumsum(oct_counts)))
fc = oct_handler.fcoords(domain,mask,cell_count, level_counts.copy())
leavesb = oct_handler.count_leaves(mask)
assert leaves == leavesb

#Now take the fcoords, call them particles and recreate the same octree
print("particle-based recreate")
oct_handler2 = ParticleOctreeContainer(dd,dle,dre)
oct_handler2.allocate_domains([nocts])
oct_handler2.n_ref = 1 #specifically make a maximum of 1 particle per oct
oct_handler2.add(fc, 1)
print("added particles")
cell_count2 = nocts*8
oct_counts2 = oct_handler.count_levels(max_level, 1, mask)
level_counts2 = np.concatenate(([0,],np.cumsum(oct_counts)))
fc2 = oct_handler.fcoords(domain,mask,cell_count, level_counts.copy())
leaves2 = oct_handler2.count_leaves(mask)
assert leaves == leaves2

print("success")

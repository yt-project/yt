fKD

****** fKD_module vars:
t1 _tree_set
t2 _tree_set
# This contains the objects that any given kD tree might wish to use.
# Not all of these are being used all the time, but they only take memory
# if they're initialized in python.
tags(:) _integer # particle ID tags
dist(:) _real # interparticle spacings
nn_tags(:,:) _integer # for all particles at once, [nth neighbor, index]
chunk_tags(:,:) _integer # for finding only a chunk of the nearest neighbors
nn_dist(:,:) _real 
pos(3,:) _real
dens(:) _real
mass(:) _real
qv(3) real
qv_many(3,:) _real
nparts integer
nn integer
nMerge integer # number of nearest neighbors used in chain merging
start integer
finish integer
tree2 _kdtree2
sort logical /.false./
rearrange logical /.true./
radius real # the unsquared radius for radius searches
radius_n integer # the number of results to return
nfound integer # number of neighbors within radius
nfound_many(:) _integer # an array of number of neighbors within radius

%%%% interval:
lower real
upper real
#real(kdkind) :: lower,upper


%%%% tree_node:
# an internal tree node
cut_dim integer
#integer :: cut_dim
# the dimension to cut
cut_val real
#real(kdkind) :: cut_val
# where to cut the dimension
cut_val_left real
cut_val_right real
#real(kdkind) :: cut_val_left, cut_val_right  
# improved cutoffs knowing the spread in child boxes.
u integer
l integer
#integer :: l, u
left _tree_node
right _tree_node
#type(tree_node), pointer :: left, right
box(:) _interval
#type(interval), pointer :: box(:) => null()
# child pointers
# Points included in this node are indexes[k] with k \in [l,u] 


%%%% kdtree2:
# Global information about the tree, one per tree
dimen integer /0/
n integer /0/
# dimensionality and total # of points
the_data(:,:) _real
#real(kdkind), pointer :: the_data(:,:) => null()
# pointer to the actual data array 
# 
#  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
#  which may be opposite of what may be conventional.
#  This is, because in Fortran, the memory layout is such that
#  the first dimension is in sequential order.  Hence, with
#  (1:d,1:N), all components of the vector will be in consecutive
#  memory locations.  The search time is dominated by the
#  evaluation of distances in the terminal nodes.  Putting all
#  vector components in consecutive memory location improves
#  memory cache locality, and hence search speed, and may enable 
#  vectorization on some processors and compilers. 
ind(:) _integer
#integer, pointer :: ind(:) => null()
# permuted index into the data, so that indexes[l..u] of some
# bucket represent the indexes of the actual points in that
# bucket.
# do we always sort output results?
sort logical /.false./
#logical       :: sort = .false.
rearrange logical /.false./
#logical       :: rearrange = .false. 
rearranged_data(:,:) _real
#real(kdkind), pointer :: rearranged_data(:,:) => null()
# if (rearrange .eqv. .true.) then rearranged_data has been
# created so that rearranged_data(:,i) = the_data(:,ind(i)),
# permitting search to use more cache-friendly rearranged_data, at
# some initial computation and storage cost.
root _tree_node
#type(tree_node), pointer :: root => null()
# root pointer of the tree


%%%% tree_set:
# This contains the objects that any given kD tree might wish to use.
# Not all of these are being used all the time, but they only take memory
# if they're initialized in python.
tags(:) _integer # particle ID tags
dist(:) _real # interparticle spacings
nn_tags(:,:) _integer # for all particles at once, [nth neighbor, index]
chunk_tags(:,:) _integer # for finding only a chunk of the nearest neighbors
nn_dist(:,:) _real 
pos(3,:) _real
dens(:) _real
mass(:) _real
qv(3) real
qv_many(3,:) _real
nparts integer
nn integer
nMerge integer # number of nearest neighbors used in chain merging
start integer
finish integer
tree2 _kdtree2
sort logical /.false./
rearrange logical /.true./
radius real # the unsquared radius for radius searches
radius_n integer # the number of results to return
nfound integer # number of neighbors within radius
nfound_many(:) _integer # an array of number of neighbors within radius

***** Subroutines:
find_nn_nearest_neighbors subroutine
create_tree(treeID:integer) subroutine
add_tree(treeID:integer) subroutine
free_tree(treeID:integer) subroutine
find_all_nn_nearest_neighbors subroutine
find_r_nearest subroutine
find_many_r_nearest(treeID:integer) subroutine
find_many_nn_nearest_neighbors subroutine
find_chunk_nearest_neighbors subroutine
chainHOP_tags_dens subroutine

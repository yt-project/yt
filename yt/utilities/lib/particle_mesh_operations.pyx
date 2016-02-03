"""
Simple integrators for the radiative transfer equation



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
cimport cython
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip

@cython.boundscheck(False)
@cython.wraparound(False)
def CICDeposit_3(np.ndarray[np.float64_t, ndim=1] posx,
                 np.ndarray[np.float64_t, ndim=1] posy,
                 np.ndarray[np.float64_t, ndim=1] posz,
                 np.ndarray[np.float64_t, ndim=1] mass,
                 np.int64_t npositions,
                 np.ndarray[np.float64_t, ndim=3] field,
                 np.ndarray[np.float64_t, ndim=1] leftEdge,
                 np.ndarray[np.int32_t, ndim=1] gridDimension,
                 np.float64_t cellSize):

    cdef int i1, j1, k1, n
    cdef np.float64_t xpos, ypos, zpos
    cdef np.float64_t fact, edge0, edge1, edge2
    cdef np.float64_t le0, le1, le2
    cdef np.float64_t dx, dy, dz, dx2, dy2, dz2

    edge0 = (<np.float64_t> gridDimension[0]) - 0.5001
    edge1 = (<np.float64_t> gridDimension[1]) - 0.5001
    edge2 = (<np.float64_t> gridDimension[2]) - 0.5001
    fact = 1.0 / cellSize

    le0 = leftEdge[0]
    le1 = leftEdge[1]
    le2 = leftEdge[2]

    for n in range(npositions):

        # Compute the position of the central cell
        xpos = (posx[n] - le0)*fact
        ypos = (posy[n] - le1)*fact
        zpos = (posz[n] - le2)*fact

        if (xpos < 0.5001) or (xpos > edge0):
            continue
        if (ypos < 0.5001) or (ypos > edge1):
            continue
        if (zpos < 0.5001) or (zpos > edge2):
            continue

        i1  = <int> (xpos + 0.5)
        j1  = <int> (ypos + 0.5)
        k1  = <int> (zpos + 0.5)

        # Compute the weights
        dx = (<np.float64_t> i1) + 0.5 - xpos
        dy = (<np.float64_t> j1) + 0.5 - ypos
        dz = (<np.float64_t> k1) + 0.5 - zpos
        dx2 =  1.0 - dx
        dy2 =  1.0 - dy
        dz2 =  1.0 - dz

        # Interpolate from field into sumfield
        field[i1-1,j1-1,k1-1] += mass[n] * dx  * dy  * dz
        field[i1  ,j1-1,k1-1] += mass[n] * dx2 * dy  * dz
        field[i1-1,j1  ,k1-1] += mass[n] * dx  * dy2 * dz
        field[i1  ,j1  ,k1-1] += mass[n] * dx2 * dy2 * dz
        field[i1-1,j1-1,k1  ] += mass[n] * dx  * dy  * dz2
        field[i1  ,j1-1,k1  ] += mass[n] * dx2 * dy  * dz2
        field[i1-1,j1  ,k1  ] += mass[n] * dx  * dy2 * dz2
        field[i1  ,j1  ,k1  ] += mass[n] * dx2 * dy2 * dz2

@cython.boundscheck(False)
@cython.wraparound(False)
def CICDeposit_2(np.ndarray[np.float64_t, ndim=1] posx,
                 np.ndarray[np.float64_t, ndim=1] posy,
                 np.ndarray[np.float64_t, ndim=1] mass,
                 np.int64_t npositions,
                 np.ndarray[np.float64_t, ndim=2] field,
                 np.ndarray[np.float64_t, ndim=1] leftEdge,
                 np.ndarray[np.int32_t, ndim=1] gridDimension,
                 np.ndarray[np.float64_t, ndim=1] cellSize):

    cdef int i1, j1, n
    cdef np.float64_t xpos, ypos
    cdef np.float64_t edge0, edge1
    cdef np.float64_t le0, le1
    cdef np.float64_t dx, dy, dx2, dy2

    edge0 = (<np.float64_t> gridDimension[0]) - 0.5001
    edge1 = (<np.float64_t> gridDimension[1]) - 0.5001

    le0 = leftEdge[0]
    le1 = leftEdge[1]

    for n in range(npositions):

        # Compute the position of the central cell
        xpos = (posx[n] - le0)/cellSize[0]
        ypos = (posy[n] - le1)/cellSize[1]

        if (xpos < 0.5001) or (xpos > edge0):
            continue
        if (ypos < 0.5001) or (ypos > edge1):
            continue

        i1  = <int> (xpos + 0.5)
        j1  = <int> (ypos + 0.5)

        # Compute the weights
        dx = (<np.float64_t> i1) + 0.5 - xpos
        dy = (<np.float64_t> j1) + 0.5 - ypos
        dx2 =  1.0 - dx
        dy2 =  1.0 - dy

        # Deposit onto field
        field[i1-1,j1-1] += mass[n] * dx  * dy
        field[i1  ,j1-1] += mass[n] * dx2 * dy
        field[i1-1,j1  ] += mass[n] * dx  * dy2
        field[i1  ,j1  ] += mass[n] * dx2 * dy2

@cython.boundscheck(False)
@cython.wraparound(False)
def NGPDeposit_2(np.ndarray[np.float64_t, ndim=1] posx,
                 np.ndarray[np.float64_t, ndim=1] posy,
                 np.ndarray[np.float64_t, ndim=1] mass,
                 np.int64_t npositions,
                 np.ndarray[np.float64_t, ndim=2] field,
                 np.ndarray[np.float64_t, ndim=1] leftEdge,
                 np.ndarray[np.int32_t, ndim=1] gridDimension,
                 np.ndarray[np.float64_t, ndim=1] cellSize):

    cdef int i1, j1, n
    cdef np.float64_t xpos, ypos
    cdef np.float64_t edge0, edge1
    cdef np.float64_t le0, le1

    edge0 = (<np.float64_t> gridDimension[0]) - 0.5001
    edge1 = (<np.float64_t> gridDimension[1]) - 0.5001

    le0 = leftEdge[0]
    le1 = leftEdge[1]

    for n in range(npositions):

        # Compute the position of the central cell
        xpos = (posx[n] - le0)/cellSize[0]
        ypos = (posy[n] - le1)/cellSize[1]

        if (xpos < 0.5001) or (xpos > edge0):
            continue
        if (ypos < 0.5001) or (ypos > edge1):
            continue

        i1  = <int> (xpos + 0.5)
        j1  = <int> (ypos + 0.5)

        # Deposit onto field
        field[i1,j1] += mass[n]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sample_field_at_positions(np.ndarray[np.float64_t, ndim=3] arr,
                              np.ndarray[np.float64_t, ndim=1] left_edge,
                              np.ndarray[np.float64_t, ndim=1] right_edge,
                              np.ndarray[np.float64_t, ndim=1] pos_x,
                              np.ndarray[np.float64_t, ndim=1] pos_y,
                              np.ndarray[np.float64_t, ndim=1] pos_z):
    cdef np.float64_t idds[3]
    cdef int dims[3]
    cdef int ind[3]
    cdef int i, npart
    npart = pos_x.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] sample
    sample = np.zeros(npart, dtype='float64')
    for i in range(3):
        dims[i] = arr.shape[i]
        idds[i] = (<np.float64_t> dims[i]) / (right_edge[i] - left_edge[i])
    for i in range(npart):
        if not ((left_edge[0] <= pos_x[i] <= right_edge[0]) and
                (left_edge[1] <= pos_y[i] <= right_edge[1]) and
                (left_edge[2] <= pos_z[i] <= right_edge[2])):
            continue
        ind[0] = <int> ((pos_x[i] - left_edge[0]) * idds[0])
        ind[1] = <int> ((pos_y[i] - left_edge[1]) * idds[1])
        ind[2] = <int> ((pos_z[i] - left_edge[2]) * idds[2])
        sample[i] = arr[ind[0], ind[1], ind[2]]
    return sample

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def CICSample_3(np.ndarray[np.float64_t, ndim=1] posx,
                np.ndarray[np.float64_t, ndim=1] posy,
                np.ndarray[np.float64_t, ndim=1] posz,
                np.ndarray[np.float64_t, ndim=1] sample,
                np.int64_t npositions,
                np.ndarray[np.float64_t, ndim=3] field,
                np.ndarray[np.float64_t, ndim=1] leftEdge,
                np.ndarray[np.int32_t, ndim=1] gridDimension,
                np.float64_t cellSize):

    cdef int i1, j1, k1, n
    cdef np.float64_t xpos, ypos, zpos
    cdef np.float64_t fact, edge0, edge1, edge2
    cdef np.float64_t le0, le1, le2
    cdef np.float64_t dx, dy, dz, dx2, dy2, dz2

    edge0 = (<np.float64_t> gridDimension[0]) - 0.5001
    edge1 = (<np.float64_t> gridDimension[1]) - 0.5001
    edge2 = (<np.float64_t> gridDimension[2]) - 0.5001
    fact = 1.0 / cellSize

    le0 = leftEdge[0]
    le1 = leftEdge[1]
    le2 = leftEdge[2]

    for n in range(npositions):

        # Compute the position of the central cell

        xpos = (posx[n]-le0)*fact
        ypos = (posy[n]-le1)*fact
        zpos = (posz[n]-le2)*fact

        if (xpos < -1 or ypos < -1 or zpos < -1 or
            xpos >= edge0+1.5001 or ypos >= edge1+1.5001 or zpos >= edge2+1.5001):
            continue

        xpos = fclip(xpos, 0.5001, edge0)
        ypos = fclip(ypos, 0.5001, edge1)
        zpos = fclip(zpos, 0.5001, edge2)

        i1  = <int> (xpos + 0.5)
        j1  = <int> (ypos + 0.5)
        k1  = <int> (zpos + 0.5)

        # Compute the weights
        dx = (<float> i1) + 0.5 - xpos
        dy = (<float> j1) + 0.5 - ypos
        dz = (<float> k1) + 0.5 - zpos
        dx2 =  1.0 - dx
        dy2 =  1.0 - dy
        dz2 =  1.0 - dz

        # Interpolate from field onto the particle
        sample[n] = (field[i1-1,j1-1,k1-1] * dx  * dy  * dz +
                     field[i1  ,j1-1,k1-1] * dx2 * dy  * dz +
                     field[i1-1,j1  ,k1-1] * dx  * dy2 * dz +
                     field[i1  ,j1  ,k1-1] * dx2 * dy2 * dz +
                     field[i1-1,j1-1,k1  ] * dx  * dy  * dz2 +
                     field[i1  ,j1-1,k1  ] * dx2 * dy  * dz2 +
                     field[i1-1,j1  ,k1  ] * dx  * dy2 * dz2 +
                     field[i1  ,j1  ,k1  ] * dx2 * dy2 * dz2)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def assign_particles_to_cells(np.ndarray[np.int32_t, ndim=1] levels, #for cells
                              np.ndarray[np.float32_t, ndim=2] left_edges, #many cells
                              np.ndarray[np.float32_t, ndim=2] right_edges,
                              np.ndarray[np.float32_t, ndim=1] pos_x, #particle
                              np.ndarray[np.float32_t, ndim=1] pos_y,
                              np.ndarray[np.float32_t, ndim=1] pos_z):
    #for every cell, assign the particles belonging to it,
    #skipping previously assigned particles
    cdef long level_max = np.max(levels)
    cdef long i,j,level
    cdef long npart = pos_x.shape[0]
    cdef long ncells = left_edges.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] assign = np.zeros(npart,dtype='int32')-1
    for level in range(level_max,0,-1):
        #start with the finest level
        for i in range(ncells):
            #go through every cell on the finest level first
            if not levels[i] == level: continue
            for j in range(npart):
                #iterate over all particles, skip if assigned
                if assign[j]>-1: continue
                if (left_edges[i,0] <= pos_x[j] <= right_edges[i,0]):
                    if (left_edges[i,1] <= pos_y[j] <= right_edges[i,1]):
                        if (left_edges[i,2] <= pos_z[j] <= right_edges[i,2]):
                            assign[j]=i
    return assign



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def assign_particles_to_cell_lists(np.ndarray[np.int32_t, ndim=1] levels, #for cells
                              np.ndarray[np.int32_t,ndim=1] assign,
                              np.int64_t level_max,
                              np.ndarray[np.float32_t, ndim=2] left_edges, #many cells
                              np.ndarray[np.float32_t, ndim=2] right_edges,
                              np.ndarray[np.float32_t, ndim=1] pos_x, #particle
                              np.ndarray[np.float32_t, ndim=1] pos_y,
                              np.ndarray[np.float32_t, ndim=1] pos_z):
    #for every cell, assign the particles belonging to it,
    #skipping previously assigned particles
    #Todo: instead of iterating every particles, could use kdtree
    cdef long i,j,level
    cdef long npart = pos_x.shape[0]
    cdef long ncells = left_edges.shape[0]
    #cdef np.ndarray[np.int32_t, ndim=1] assign
    #assign = np.zeros(npart,dtype='int32')-1
    index_lists = []
    for level in range(level_max,-1,-1):
        #start with the finest level
        for i in range(ncells):
            #go through every cell on the finest level first
            if not levels[i] == level: continue
            index_list = []
            for j in range(npart):
                #iterate over all particles, skip if assigned
                if assign[j]>-1: continue
                if (left_edges[i,0] <= pos_x[j] <= right_edges[i,0]):
                    if (left_edges[i,1] <= pos_y[j] <= right_edges[i,1]):
                        if (left_edges[i,2] <= pos_z[j] <= right_edges[i,2]):
                            assign[j]=i
                            index_list += j,
            index_lists += index_list,
    return assign,index_lists


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def recursive_particle_assignment(grids, grid,
                                  np.ndarray[np.float32_t, ndim=2] left_edges, #many cells
                                  np.ndarray[np.float32_t, ndim=2] right_edges,
                                  np.ndarray[np.float32_t, ndim=1] pos_x, #particle
                                  np.ndarray[np.float32_t, ndim=1] pos_y,
                                  np.ndarray[np.float32_t, ndim=1] pos_z):
    #start on level zero, grid particles onto every mesh
    #every particle we are fed, we can assume it exists on our grid
    #must fill in the grid_particle_count array
    #and particle_indices for every grid
    cdef long i, j
    cdef long npart = pos_x.shape[0]
    cdef np.ndarray[np.int32_t, ndim=1] assigned       = np.zeros(npart,dtype='int32')
    cdef np.ndarray[np.int32_t, ndim=1] never_assigned = np.ones(npart,dtype='int32')
    for i in np.unique(grid.child_index_mask):
        if i== -1: continue
        #assigned to this subgrid
        assigned = np.zeros(npart,dtype='int32')
        for j in range(npart):
            if (left_edges[i,0] <= pos_x[j] <= right_edges[i,0]):
                if (left_edges[i,1] <= pos_y[j] <= right_edges[i,1]):
                    if (left_edges[i,2] <= pos_z[j] <= right_edges[i,2]):
                       assigned[j]=1
                       never_assigned[j]=0
        if np.sum(assigned)>0:
            recursive_particle_assignment(grids,grid,left_edges,right_edges,
                                           pos_x[assigned],pos_y[assigned],pos_z[assigned])
    #now we have assigned particles to other subgrids, we are left with particles on our grid






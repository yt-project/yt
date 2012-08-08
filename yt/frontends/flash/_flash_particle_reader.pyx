import numpy as np
cimport numpy as np
cimport cython
import h5py

cdef particles_validator_region(np.ndarray[np.float64_t, ndim=1] x,
                                np.ndarray[np.float64_t, ndim=1] y,
                                np.ndarray[np.float64_t, ndim=1] z,
                                np.ndarray[np.float64_t, ndim=1] left_edge,
                                np.ndarray[np.float64_t, ndim=1] right_edge,
                                np.int32_t periodic,
                                np.ndarray[np.float64_t, ndim=1] DLE,
                                np.ndarray[np.float64_t, ndim=1] DRE) :

    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxs
    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxx
    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxy
    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxz
    
    cdef np.ndarray[np.float64_t, ndim=1] xx
    cdef np.ndarray[np.float64_t, ndim=1] yy
    cdef np.ndarray[np.float64_t, ndim=1] zz

    cdef np.ndarray[np.float64_t, ndim=1] DW

    idxs = np.zeros(x.shape[0], 'bool')
    idxx = np.zeros(x.shape[0], 'bool')
    idxy = np.zeros(x.shape[0], 'bool')
    idxz = np.zeros(x.shape[0], 'bool')

    xx = np.zeros(x.shape[0], 'float64')
    yy = np.zeros(x.shape[0], 'float64')
    zz = np.zeros(x.shape[0], 'float64')

    DW = np.zeros(3, 'float64')

    xx = x
    yy = y
    zz = z

    if periodic == 1 : 

        DW = DRE - DLE
        xx[x < left_edge[0]] = x + DW[0]
        xx[x > right_edge[0]] = x - DW[0]
        yy[y < left_edge[1]] = y + DW[1]
        yy[y > right_edge[1]] = y - DW[1]
        zz[z < left_edge[2]] = z + DW[2]
        zz[z > right_edge[2]] = z - DW[2]

    idxx = np.logical_and(xx >= left_edge[0], xx <= right_edge[0])
    idxy = np.logical_and(yy >= left_edge[1], yy <= right_edge[1])
    idxz = np.logical_and(zz >= left_edge[2], zz <= right_edge[2])

    idxs = np.logical_and(idxx, idxy)
    idxs = np.logical_and(idxz, idxs)

    return idxs

cdef particles_validator_sphere(np.ndarray[np.float64_t, ndim=1] x,
                                np.ndarray[np.float64_t, ndim=1] y, 
                                np.ndarray[np.float64_t, ndim=1] z,
                                np.ndarray[np.float64_t, ndim=1] center,
                                np.float64_t radius,
                                np.int32_t periodic,
                                np.ndarray[np.float64_t, ndim=1] DLE,
                                np.ndarray[np.float64_t, ndim=1] DRE) :

    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxs

    cdef np.ndarray[np.float64_t, ndim=1] r
    cdef np.ndarray[np.float64_t, ndim=1] xx
    cdef np.ndarray[np.float64_t, ndim=1] yy
    cdef np.ndarray[np.float64_t, ndim=1] zz

    cdef np.ndarray[np.float64_t, ndim=1] DW

    idxs = np.zeros(x.shape[0], 'bool')
    
    r = np.zeros(x.shape[0], 'float64')
    xx = np.zeros(x.shape[0], 'float64')
    yy = np.zeros(x.shape[0], 'float64')
    zz = np.zeros(x.shape[0], 'float64')

    DW = np.zeros(3, 'float64')
    
    xx = np.abs(x-center[0])
    yy = np.abs(y-center[1])
    zz = np.abs(z-center[2])

    if periodic == 1 : 

        DW = DRE - DLE

        xx = np.minimum(xx,DW[0]-xx)
        yy = np.minimum(yy,DW[1]-yy)
        zz = np.minimum(zz,DW[2]-zz)

    r = np.sqrt(xx*xx+yy*yy+zz*zz)

    idxs = np.array(r <= radius)
    
    return idxs

cdef particles_validator_disk(np.ndarray[np.float64_t, ndim=1] x,
                              np.ndarray[np.float64_t, ndim=1] y,
                              np.ndarray[np.float64_t, ndim=1] z,
                              np.ndarray[np.float64_t, ndim=1] center,
                              np.ndarray[np.float64_t, ndim=1] normal,
                              np.float64_t radius, np.float64_t height) :

    cdef np.float64_t d

    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxs

    cdef np.ndarray[np.float64_t, ndim=1] ph
    cdef np.ndarray[np.float64_t, ndim=1] pd2
    cdef np.ndarray[np.float64_t, ndim=1] pr

    idxs = np.zeros(x.shape[0], 'bool')
    
    ph = np.zeros(x.shape[0], 'float64')
    pd2 = np.zeros(x.shape[0], 'float64')
    pr = np.zeros(x.shape[0], 'float64')
    
    d = -np.dot(normal*center)

    ph = np.abs(x*normal[0] + y*normal[1] + z*normal[2] + d)
    pd2 = (x-center[0])**2+(y-center[1])**2+(z-center[2])**2

    pr = np.sqrt(pd2-ph*ph)

    idxs = np.logical_and(pr <= radius, ph <= height)
    
    return idxs

@cython.boundscheck(False)
@cython.wraparound(False)
def read_particles(file_id, int x_index, int y_index, int z_index,
                   int num_fields, int rtype, args,
                   np.ndarray[np.int32_t, ndim=1] field_indices) :

    cdef np.ndarray[np.uint8_t, cast=True, ndim=1] idxs
    cdef int i
    cdef int num_particles
    cdef np.int32_t periodic
    cdef np.ndarray[np.float64_t, ndim=1] left_edge
    cdef np.ndarray[np.float64_t, ndim=1] right_edge
    cdef np.ndarray[np.float64_t, ndim=1] DLE
    cdef np.ndarray[np.float64_t, ndim=1] DRE
    cdef np.float64_t radius
    cdef np.float64_t height
    cdef np.ndarray[np.float64_t, ndim=1] normal
    cdef np.ndarray[np.float64_t, ndim=1] center
    cdef np.ndarray[np.float64_t, ndim=1] particle_field
    cdef np.ndarray[np.float64_t, ndim=1] posx
    cdef np.ndarray[np.float64_t, ndim=1] posy
    cdef np.ndarray[np.float64_t, ndim=1] posz

    left_edge = np.zeros(3, 'float64')
    right_edge = np.zeros(3, 'float64')
    DLE = np.zeros(3, 'float64')
    DRE = np.zeros(3, 'float64')
    normal = np.zeros(3, 'float64')
    center = np.zeros(3, 'float64')

    dataset = h5py.h5d.open(file_id, "tracer particles")
    dataspace = dataset.get_space()
    rank = dataspace.get_simple_extent_dims()
    memspace = h5py.h5s.create_simple((rank[0],))

    num_particles = rank[0]
    count = (num_particles,1)

    posx = np.zeros(num_particles, 'float64')
    posy = np.zeros(num_particles, 'float64')
    posz = np.zeros(num_particles, 'float64')

    start = (0,x_index)
    dataspace.select_hyperslab(start,count)
    dataset.read(memspace, dataspace, posx)

    start = (0,y_index)
    dataspace.select_hyperslab(start,count)
    dataset.read(memspace, dataspace, posy)

    start = (0,z_index)
    dataspace.select_hyperslab(start,count)
    dataset.read(memspace, dataspace, posz)
    
    idxs = np.zeros(num_particles, 'bool')

    particle_field = np.zeros(num_particles, 'float64')
    
    if rtype == 0 :
        left_edge = args[0]
        right_edge = args[1]
        periodic = args[2]
        DLE = args[3]
        DRE = args[4]
        idxs = particles_validator_region(posx,posy,posz,
                                          left_edge,right_edge,
                                          periodic,DLE,DRE)
    elif rtype == 1:
        center = args[0]
        radius = args[1]
        periodic = args[2]
        DLE = args[3]
        DRE = args[4]
        idxs = particles_validator_sphere(posx,posy,posz,
                                          center,radius,
                                          periodic,DLE,DRE)
    elif rtype == 2:
        center = args[0]
        normal = args[1]
        radius = args[2]
        height = args[3]
        idxs = particles_validator_disk(posx,posy,posz,
                                        center,normal,
                                        radius,height)

    _particles = []

    for i in range(num_fields) :

        start = (0,field_indices[i])
        dataspace.select_hyperslab(start,count)
        dataset.read(memspace, dataspace, particle_field)
        _particles.append(particle_field[idxs])
        
    return _particles
    

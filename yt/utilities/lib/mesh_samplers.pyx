"""
This file contains coordinate mappings between physical coordinates and those
defined on unit elements, as well as functions that do the corresponding intracell
interpolation.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
from pyembree.rtcore cimport Vec3f, Triangle, Vertex
from yt.utilities.lib.mesh_construction cimport MeshDataContainer
cimport numpy as np
cimport cython
from libc.math cimport fabs, fmax


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double determinant_3x3(double* col0, 
                            double* col1, 
                            double* col2) nogil:
    return col0[0]*col1[1]*col2[2] - col0[0]*col1[2]*col2[1] - \
           col0[1]*col1[0]*col2[2] + col0[1]*col1[2]*col2[0] + \
           col0[2]*col1[0]*col2[1] - col0[2]*col1[1]*col2[0]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double maxnorm(double* f) nogil:
    cdef double err
    cdef int i
    err = fabs(f[0])
    for i in range(1, 2):
        err = fmax(err, fabs(f[i])) 
    return err


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void get_hit_position(double* position,
                           void* userPtr,
                           rtcr.RTCRay& ray) nogil:
    cdef int primID, i
    cdef double[3][3] vertex_positions
    cdef Triangle tri
    cdef MeshDataContainer* data

    primID = ray.primID
    data = <MeshDataContainer*> userPtr
    tri = data.indices[primID]

    vertex_positions[0][0] = data.vertices[tri.v0].x
    vertex_positions[0][1] = data.vertices[tri.v0].y
    vertex_positions[0][2] = data.vertices[tri.v0].z

    vertex_positions[1][0] = data.vertices[tri.v1].x
    vertex_positions[1][1] = data.vertices[tri.v1].y
    vertex_positions[1][2] = data.vertices[tri.v1].z

    vertex_positions[2][0] = data.vertices[tri.v2].x
    vertex_positions[2][1] = data.vertices[tri.v2].y
    vertex_positions[2][2] = data.vertices[tri.v2].z

    for i in range(3):
        position[i] = vertex_positions[0][i]*(1.0 - ray.u - ray.v) + \
                      vertex_positions[1][i]*ray.u + \
                      vertex_positions[2][i]*ray.v


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void linear_hex_f(double* f,
                              double* x, 
                              double* vertices, 
                              double* phys_x) nogil:
    
    cdef int i
    cdef double rm, rp, sm, sp, tm, tp
    
    rm = 1.0 - x[0]
    rp = 1.0 + x[0]
    sm = 1.0 - x[1]
    sp = 1.0 + x[1]
    tm = 1.0 - x[2]
    tp = 1.0 + x[2]
    
    for i in range(3):
        f[i] = vertices[0 + i]*rm*sm*tm \
             + vertices[3 + i]*rp*sm*tm \
             + vertices[6 + i]*rp*sp*tm \
             + vertices[9 + i]*rm*sp*tm \
             + vertices[12 + i]*rm*sm*tp \
             + vertices[15 + i]*rp*sm*tp \
             + vertices[18 + i]*rp*sp*tp \
             + vertices[21 + i]*rm*sp*tp \
             - 8.0*phys_x[i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void linear_hex_J(double* r,
                              double* s,
                              double* t,
                              double* x, 
                              double* v, 
                              double* phys_x) nogil:
    
    cdef int i
    cdef double rm, rp, sm, sp, tm, tp
    
    rm = 1.0 - x[0]
    rp = 1.0 + x[0]
    sm = 1.0 - x[1]
    sp = 1.0 + x[1]
    tm = 1.0 - x[2]
    tp = 1.0 + x[2]
    
    for i in range(3):
        r[i] = -sm*tm*v[0 + i]  + sm*tm*v[3 + i]  + \
                sp*tm*v[6 + i]  - sp*tm*v[9 + i]  - \
                sm*tp*v[12 + i] + sm*tp*v[15 + i] + \
                sp*tp*v[18 + i] - sp*tp*v[21 + i]
        s[i] = -rm*tm*v[0 + i]  - rp*tm*v[3 + i]  + \
                rp*tm*v[6 + i]  + rm*tm*v[9 + i]  - \
                rm*tp*v[12 + i] - rp*tp*v[15 + i] + \
                rp*tp*v[18 + i] + rm*tp*v[21 + i]
        t[i] = -rm*sm*v[0 + i]  - rp*sm*v[3 + i]  - \
                rp*sp*v[6 + i]  - rm*sp*v[9 + i]  + \
                rm*sm*v[12 + i] + rp*sm*v[15 + i] + \
                rp*sp*v[18 + i] + rm*sp*v[21 + i]
                
                
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double sample_hex_at_unit_point(double* coord, double* vals) nogil:
    cdef double F, rm, rp, sm, sp, tm, tp
    
    rm = 1.0 - coord[0]
    rp = 1.0 + coord[0]
    sm = 1.0 - coord[1]
    sp = 1.0 + coord[1]
    tm = 1.0 - coord[2]
    tp = 1.0 + coord[2]
    
    F = vals[0]*rm*sm*tm + vals[1]*rp*sm*tm + vals[2]*rp*sp*tm + vals[3]*rm*sp*tm + \
        vals[4]*rm*sm*tp + vals[5]*rp*sm*tp + vals[6]*rp*sp*tp + vals[7]*rm*sp*tp
    return 0.125*F
                

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int hex_check_inside(double* mapped_coord) nogil:    
    if (fabs(mapped_coord[0]) - 1.0 > 1.0e-8 or
        fabs(mapped_coord[1]) - 1.0 > 1.0e-8 or 
        fabs(mapped_coord[2]) - 1.0 > 1.0e-8):
        return 0
    return 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double sample_hex_at_real_point(double* vertices,
                                     double* field_values,
                                     double* physical_x) nogil:
    
    cdef int i
    cdef double d, val
    cdef double[3] f
    cdef double[3] r
    cdef double[3] s
    cdef double[3] t
    cdef double[3] x
    cdef double tolerance = 1.0e-9
    cdef int iterations = 0
    cdef double err
   
    # initial guess
    for i in range(3):
        x[i] = 0.0
    
    # initial error norm
    linear_hex_f(f, x, vertices, physical_x)
    err = maxnorm(f)  
   
    # begin Newton iteration
    while (err > tolerance and iterations < 100):
        linear_hex_J(r, s, t, x, vertices, physical_x)
        d = determinant_3x3(r, s, t)
        x[0] = x[0] - (determinant_3x3(f, s, t)/d)
        x[1] = x[1] - (determinant_3x3(r, f, t)/d)
        x[2] = x[2] - (determinant_3x3(r, s, f)/d)
        linear_hex_f(f, x, vertices, physical_x)        
        err = maxnorm(f)
        iterations += 1

    val = sample_hex_at_unit_point(x, field_values)
    return val

    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void sample_hex(void* userPtr,
                     rtcr.RTCRay& ray) nogil:
    cdef int ray_id, elem_id, i
    cdef double val
    cdef double[8] field_data
    cdef int[8] element_indices
    cdef double[24] vertices
    cdef double[3] position
    cdef MeshDataContainer* data

    data = <MeshDataContainer*> userPtr
    ray_id = ray.primID
    if ray_id == -1:
        return

    elem_id = ray_id / data.tpe

    get_hit_position(position, userPtr, ray)
    
    for i in range(8):
        element_indices[i] = data.element_indices[elem_id*8+i]
        field_data[i]      = data.field_data[elem_id*8+i]

    for i in range(8):
        vertices[i*3]     = data.vertices[element_indices[i]].x
        vertices[i*3 + 1] = data.vertices[element_indices[i]].y
        vertices[i*3 + 2] = data.vertices[element_indices[i]].z    

    val = sample_hex_at_real_point(vertices, field_data, position)
    ray.time = val


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void hex_real_to_mapped(double* mapped_x,
                             double* vertices,
                             double* physical_x) nogil:
    
    cdef int i
    cdef double d, val
    cdef double[3] f
    cdef double[3] r
    cdef double[3] s
    cdef double[3] t
    cdef double[3] x
    cdef double tolerance = 1.0e-9
    cdef int iterations = 0
    cdef double err
   
    # initial guess
    for i in range(3):
        x[i] = 0.0
    
    # initial error norm
    linear_hex_f(f, x, vertices, physical_x)
    err = maxnorm(f)  
   
    # begin Newton iteration
    while (err > tolerance and iterations < 10):
        linear_hex_J(r, s, t, x, vertices, physical_x)
        d = determinant_3x3(r, s, t)
        x[0] = x[0] - (determinant_3x3(f, s, t)/d)
        x[1] = x[1] - (determinant_3x3(r, f, t)/d)
        x[2] = x[2] - (determinant_3x3(r, s, f)/d)
        linear_hex_f(f, x, vertices, physical_x)        
        err = maxnorm(f)
        iterations += 1

    for i in range(3):
        mapped_x[i] = x[i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def test_hex_sampler(np.ndarray[np.float64_t, ndim=2] vertices,
                     np.ndarray[np.float64_t, ndim=1] field_values,
                     np.ndarray[np.float64_t, ndim=1] physical_x):

    cdef double val

    val = sample_hex_at_real_point(<double*> vertices.data,
                                   <double*> field_values.data,
                                   <double*> physical_x.data)
    return val

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double sample_tetra_at_unit_point(double* coord, double* vals) nogil:
    return vals[0]*coord[0] + vals[1]*coord[1] + vals[2]*coord[2] + vals[3]*coord[3]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int tetra_check_inside(double* mapped_coord) nogil:    
    cdef int i
    for i in range(4):
        if (mapped_coord[i] < -1.0e-8 or
            mapped_coord[i] - 1.0 > 1.0e-8):
            return 0
    return 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void tetra_real_to_mapped(double* mapped_coord,
                               double* vertices,
                               double* physical_coord) nogil:
    cdef int i
    cdef double d
    cdef double[3] bvec
    cdef double[3] col0
    cdef double[3] col1
    cdef double[3] col2
    
    for i in range(3):
        bvec[i] = physical_coord[i]   - vertices[9 + i]
        col0[i] = vertices[0 + i]     - vertices[9 + i]
        col1[i] = vertices[3 + i]     - vertices[9 + i]
        col2[i] = vertices[6 + i]     - vertices[9 + i]
        
    d = determinant_3x3(col0, col1, col2)
    mapped_coord[0] = determinant_3x3(bvec, col1, col2)/d
    mapped_coord[1] = determinant_3x3(col0, bvec, col2)/d
    mapped_coord[2] = determinant_3x3(col0, col1, bvec)/d
    mapped_coord[3] = 1.0 - mapped_coord[0] - mapped_coord[1] - mapped_coord[2]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double sample_tetra_at_real_point(double* vertices,
                                       double* field_values,
                                       double* physical_x) nogil:
    cdef double val
    cdef double mapped_coord[4]

    tetra_real_to_mapped(mapped_coord, 
                         vertices,
                         physical_x)    
        
    val = sample_tetra_at_unit_point(mapped_coord, field_values)
    return val


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void sample_tetra(void* userPtr,
                       rtcr.RTCRay& ray) nogil:

    cdef int ray_id, elem_id, i
    cdef double val
    cdef double[4] field_data
    cdef int[4] element_indices
    cdef double[12] vertices
    cdef double[3] position
    cdef MeshDataContainer* data

    data = <MeshDataContainer*> userPtr
    ray_id = ray.primID
    if ray_id == -1:
        return

    get_hit_position(position, userPtr, ray)
    
    elem_id = ray_id / 4
    for i in range(4):
        element_indices[i] = data.element_indices[elem_id*4+i]
        field_data[i] = data.field_data[elem_id*4+i]
        vertices[i*3] = data.vertices[element_indices[i]].x
        vertices[i*3 + 1] = data.vertices[element_indices[i]].y
        vertices[i*3 + 2] = data.vertices[element_indices[i]].z    

    val = sample_tetra_at_real_point(vertices, field_data, position)
    ray.time = val

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def test_tetra_sampler(np.ndarray[np.float64_t, ndim=2] vertices,
                       np.ndarray[np.float64_t, ndim=1] field_values,
                       np.ndarray[np.float64_t, ndim=1] physical_x):

    cdef double val
    cdef double[4] mapped_coord
    tetra_real_to_mapped(mapped_coord, 
                         <double*> vertices.data,
                         <double*> physical_x.data)

    val = sample_tetra_at_unit_point(mapped_coord, 
                                     <double*> field_values.data)
    return val

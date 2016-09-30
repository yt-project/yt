"""
This file contains functions that sample a surface mesh at the point hit by
a ray. These can be used with pyembree in the form of "filter feedback functions."

Note - this file is only used for the Embree-accelerated ray-tracer.

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
from yt.utilities.lib.mesh_construction cimport \
    MeshDataContainer, \
    Patch
from yt.utilities.lib.primitives cimport patchSurfaceFunc
from yt.utilities.lib.element_mappings cimport \
    ElementSampler, \
    P1Sampler3D, \
    Q1Sampler3D, \
    S2Sampler3D, \
    W1Sampler3D
cimport numpy as np
cimport cython
from libc.math cimport fabs, fmax

cdef ElementSampler Q1Sampler = Q1Sampler3D()
cdef ElementSampler P1Sampler = P1Sampler3D()
cdef ElementSampler S2Sampler = S2Sampler3D()
cdef ElementSampler W1Sampler = W1Sampler3D()


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

    # ray_id records the id number of the hit according to
    # embree, in which the primitives are triangles. Here,
    # we convert this to the element id by dividing by the
    # number of triangles per element.
    elem_id = ray_id / data.tpe

    get_hit_position(position, userPtr, ray)
    
    for i in range(8):
        element_indices[i] = data.element_indices[elem_id*8+i]

    for i in range(data.fpe):
        field_data[i] = data.field_data[elem_id*data.fpe+i]

    for i in range(8):
        vertices[i*3]     = data.vertices[element_indices[i]].x
        vertices[i*3 + 1] = data.vertices[element_indices[i]].y
        vertices[i*3 + 2] = data.vertices[element_indices[i]].z    

    # we use ray.time to pass the value of the field
    cdef double mapped_coord[3]
    Q1Sampler.map_real_to_unit(mapped_coord, vertices, position)
    if data.fpe == 1:
        val = field_data[0]
    else:
        val = Q1Sampler.sample_at_unit_point(mapped_coord, field_data)
    ray.time = val

    # we use ray.instID to pass back whether the ray is near the
    # element boundary or not (used to annotate mesh lines)
    ray.instID = Q1Sampler.check_mesh_lines(mapped_coord)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void sample_wedge(void* userPtr,
                       rtcr.RTCRay& ray) nogil:
    cdef int ray_id, elem_id, i
    cdef double val
    cdef double[6] field_data
    cdef int[6] element_indices
    cdef double[18] vertices
    cdef double[3] position
    cdef MeshDataContainer* data

    data = <MeshDataContainer*> userPtr
    ray_id = ray.primID
    if ray_id == -1:
        return

    # ray_id records the id number of the hit according to
    # embree, in which the primitives are triangles. Here,
    # we convert this to the element id by dividing by the
    # number of triangles per element.
    elem_id = ray_id / data.tpe

    get_hit_position(position, userPtr, ray)
    
    for i in range(6):
        element_indices[i] = data.element_indices[elem_id*6+i]

    for i in range(data.fpe):
        field_data[i] = data.field_data[elem_id*data.fpe+i]

    for i in range(6):
        vertices[i*3]     = data.vertices[element_indices[i]].x
        vertices[i*3 + 1] = data.vertices[element_indices[i]].y
        vertices[i*3 + 2] = data.vertices[element_indices[i]].z    

    # we use ray.time to pass the value of the field
    cdef double mapped_coord[3]
    W1Sampler.map_real_to_unit(mapped_coord, vertices, position)
    if data.fpe == 1:
        val = field_data[0]
    else:
        val = W1Sampler.sample_at_unit_point(mapped_coord, field_data)
    ray.time = val
    ray.instID = W1Sampler.check_mesh_lines(mapped_coord)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef void sample_hex20(void* userPtr,
                       rtcr.RTCRay& ray) nogil:
    cdef int ray_id, elem_id, i
    cdef double val
    cdef double[3] position
    cdef float[3] pos
    cdef Patch* data

    data = <Patch*> userPtr
    ray_id = ray.primID
    if ray_id == -1:
        return
    cdef Patch patch = data[ray_id]

    # ray_id records the id number of the hit according to
    # embree, in which the primitives are patches. Here,
    # we convert this to the element id by dividing by the
    # number of patches per element.
    elem_id = ray_id / 6

    # fills "position" with the physical position of the hit
    patchSurfaceFunc(data[ray_id].v, ray.u, ray.v, pos)
    for i in range(3):
        position[i] = <double> pos[i]
 
    # we use ray.time to pass the value of the field
    cdef double mapped_coord[3]
    S2Sampler.map_real_to_unit(mapped_coord, patch.vertices, position)
    val = S2Sampler.sample_at_unit_point(mapped_coord, patch.field_data)
    ray.time = val
    ray.instID = S2Sampler.check_mesh_lines(mapped_coord)
    

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

    # ray_id records the id number of the hit according to
    # embree, in which the primitives are triangles. Here,
    # we convert this to the element id by dividing by the
    # number of triangles per element.    
    elem_id = ray_id / data.tpe

    for i in range(4):
        element_indices[i] = data.element_indices[elem_id*4+i]
        vertices[i*3] = data.vertices[element_indices[i]].x
        vertices[i*3 + 1] = data.vertices[element_indices[i]].y
        vertices[i*3 + 2] = data.vertices[element_indices[i]].z    

    for i in range(data.fpe):
        field_data[i] = data.field_data[elem_id*data.fpe+i]

    # we use ray.time to pass the value of the field
    cdef double mapped_coord[4]
    P1Sampler.map_real_to_unit(mapped_coord, vertices, position)
    if data.fpe == 1:
        val = field_data[0]
    else:
        val = P1Sampler.sample_at_unit_point(mapped_coord, field_data)
    ray.time = val
    ray.instID = P1Sampler.check_mesh_lines(mapped_coord)

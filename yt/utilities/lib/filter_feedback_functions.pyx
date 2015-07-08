cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
from pyembree.rtcore cimport Vec3f, Triangle, Vertex
from yt.utilities.lib.mesh_construction cimport UserData
from yt.utilities.lib.element_mappings import Q1Sampler3D
cimport numpy as np
cimport cython


cdef double get_value_trilinear(void* userPtr,
                                rtcr.RTCRay& ray):
    cdef int ray_id
    cdef double u, v, val
    cdef double d0, d1, d2

    data = <UserData*> userPtr
    ray_id = ray.primID

    u = ray.u
    v = ray.v

    d0 = data.field_data[ray_id].x
    d1 = data.field_data[ray_id].y
    d2 = data.field_data[ray_id].z

    return d0*(1.0 - u - v) + d1*u + d2*v

    
cdef double sample_surface_hex(void* userPtr,
                               rtcr.RTCRay& ray):
    cdef int ray_id, elem_id, i
    cdef double u, v, val
    cdef double d0, d1, d2
    cdef double[:] field_data
    cdef long[:] element_indices
    cdef double[8][3] vertices
    cdef double[:] position
    cdef double result
    cdef UserData* data

    data = <UserData*> userPtr
    ray_id = ray.primID
    elem_id = ray_id / data.tpe

    position = get_hit_position(userPtr, ray)
    element_indices = data.element_indices[elem_id]
    field_data = data.field_data[elem_id]

    for i in range(8):
        vertices[i][0] = data.vertices[element_indices[i]].x
        vertices[i][1] = data.vertices[element_indices[i]].y
        vertices[i][2] = data.vertices[element_indices[i]].z    

    sampler = Q1Sampler3D()
    result = sampler.sample_at_real_point(position, vertices, field_data)

    return result


cdef double[:] get_hit_position(void* userPtr,
                                rtcr.RTCRay& ray):
    cdef int primID, elemID, i
    cdef double[3] position
    cdef double[3][3] vertex_positions
    cdef Triangle tri
    cdef UserData* data

    primID = ray.primID
    data = <UserData*> userPtr
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
        position[i] = vertex_positions[0][i]*ray.u + \
                      vertex_positions[1][i]*ray.v + \
                      vertex_positions[2][i]*(1.0 - ray.u - ray.v)

    return position
    

cdef void maximum_intensity(void* userPtr, 
                            rtcr.RTCRay& ray):

    cdef double val = get_value_trilinear(userPtr, ray)
    ray.time = max(ray.time, val)
    ray.geomID = -1  # reject hit


cdef void sample_surface(void* userPtr, 
                         rtcr.RTCRay& ray):

    cdef double val = get_value_trilinear(userPtr, ray)
    ray.time = val

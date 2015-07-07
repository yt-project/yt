cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
from pyembree.rtcore cimport Vec3f
from yt.utilities.lib.mesh_construction cimport UserData
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


cdef double get_value_triangle(void* userPtr,
                               rtcr.RTCRay& ray):
    cdef int ray_id, elem_id
    cdef double u, v, val
    cdef double d0, d1, d2
    cdef Vec3f* field_data
    cdef long[:, :] element_indices


    data = <UserData*> userPtr
    field_data = data.field_data
    element_indices = data.element_indices

    ray_id = ray.primID
    elem_id = ray_id / data.tpe
    u = ray.u
    v = ray.v

    d0 = field_data[ray_id].x
    d1 = field_data[ray_id].y
    d2 = field_data[ray_id].z

    return d0*(1.0 - u - v) + d1*u + d2*v


cdef void maximum_intensity(void* userPtr, 
                            rtcr.RTCRay& ray):

    cdef double val = get_value_trilinear(userPtr, ray)
    ray.time = max(ray.time, val)
    ray.geomID = -1  # reject hit


cdef void sample_surface(void* userPtr, 
                         rtcr.RTCRay& ray):

    cdef double val = get_value_trilinear(userPtr, ray)
    ray.time = val

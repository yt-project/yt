cimport pyembree.rtcore
cimport pyembree.rtcore_ray
cimport pyembree.rtcore_scene as rtcs


cdef class YTEmbreeScene:
    cdef rtcs.RTCScene scene_i

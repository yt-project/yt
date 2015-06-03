cimport pyembree.rtcore
cimport pyembree.rtcore_scene as rtcs
cimport pyembree.rtcore_ray

cdef class EmbreeVolume:
    cdef rtcs.RTCScene scene_i


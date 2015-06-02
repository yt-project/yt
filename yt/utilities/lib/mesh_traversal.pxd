cimport cython
cimport numpy as np
import numpy as np
cimport pyembree.rtcore_scene as rtcs

cdef class EmbreeVolume:
    cdef rtcs.RTCScene scene_i

"""
This file contains the MeshSampler classes, which handles casting rays at a
mesh source using either pyembree or the cython ray caster.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport cython
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_scene as rtcs
from .image_samplers cimport \
    ImageSampler
from .lenses cimport \
    ImageContainer
from cython.parallel import prange, parallel, threadid
from yt.visualization.image_writer import apply_colormap
from yt.utilities.lib.bounding_volume_hierarchy cimport BVH, Ray

rtc.rtcInit(NULL)
rtc.rtcSetErrorFunction(error_printer)

cdef void error_printer(const rtc.RTCError code, const char *_str):
    print "ERROR CAUGHT IN EMBREE"
    rtc.print_error(code)
    print "ERROR MESSAGE:", _str

cdef class YTEmbreeScene:

    def __init__(self):
        self.scene_i = rtcs.rtcNewScene(rtcs.RTC_SCENE_STATIC, rtcs.RTC_INTERSECT1)

    def __dealloc__(self):
        rtcs.rtcDeleteScene(self.scene_i)


cdef class EmbreeMeshSampler(ImageSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self, 
                 YTEmbreeScene scene,
                 int num_threads = 0):
        '''

        This function is supposed to cast the rays and return the
        image.

        '''

        rtcs.rtcCommit(scene.scene_i)
        cdef int vi, vj, i, j
        cdef ImageContainer *im = self.image
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.int64_t nx, ny, size
        cdef np.float64_t width[3]
        for i in range(3):
            width[i] = self.width[i]
        nx = im.nv[0]
        ny = im.nv[1]
        size = nx * ny
        cdef rtcr.RTCRay ray
        v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        for j in range(size):
            vj = j % ny
            vi = (j - vj) / ny
            vj = vj
            self.vector_function(im, vi, vj, width, v_dir, v_pos)
            for i in range(3):
                ray.org[i] = v_pos[i]
                ray.dir[i] = v_dir[i]
            ray.tnear = 0.0
            ray.tfar = 1e37
            ray.geomID = rtcg.RTC_INVALID_GEOMETRY_ID
            ray.primID = rtcg.RTC_INVALID_GEOMETRY_ID
            ray.instID = rtcg.RTC_INVALID_GEOMETRY_ID
            ray.mask = -1
            ray.time = 0
            ray.Ng[0] = 1e37  # we use this to track the hit distance
            rtcs.rtcIntersect(scene.scene_i, ray)
            im.image[vi, vj, 0] = ray.time
            im.image_used[vi, vj] = ray.primID
            im.mesh_lines[vi, vj] = ray.instID
            im.zbuffer[vi, vj] = ray.tfar
        free(v_pos)
        free(v_dir)

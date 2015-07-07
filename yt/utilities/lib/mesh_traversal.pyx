cimport cython
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
cimport pyembree.rtcore as rtc
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_scene as rtcs
from grid_traversal cimport ImageSampler, \
    ImageContainer
from cython.parallel import prange, parallel, threadid

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

cdef class MeshSampler(ImageSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self, 
                 YTEmbreeScene scene,
                 mesh,
                 int num_threads = 0):
        '''

        This function is supposed to cast the rays and return the
        image.

        '''

        rtcs.rtcCommit(scene.scene_i)
        cdef int vi, vj, i, j, ni, nj, nn
        cdef np.int64_t offset
        cdef ImageContainer *im = self.image
        cdef np.int64_t elemID
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.int64_t nx, ny, size
        cdef np.float64_t px, py
        cdef np.float64_t width[3]
        for i in range(3):
            width[i] = self.width[i]
        cdef np.ndarray[np.float64_t, ndim=1] data
        nx = im.nv[0]
        ny = im.nv[1]
        size = nx * ny
        data = np.empty(size, dtype="float64")
        cdef rtcr.RTCRay ray
        if im.vd_strides[0] == -1:
            v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            for j in range(size):
                vj = j % ny
                vi = (j - vj) / ny
                vj = vj
                # Dynamically calculate the position
                px = width[0] * (<np.float64_t>vi)/(<np.float64_t>im.nv[0]-1) - width[0]/2.0
                py = width[1] * (<np.float64_t>vj)/(<np.float64_t>im.nv[1]-1) - width[1]/2.0
                v_pos[0] = im.vp_pos[0]*px + im.vp_pos[3]*py + im.vp_pos[9]
                v_pos[1] = im.vp_pos[1]*px + im.vp_pos[4]*py + im.vp_pos[10]
                v_pos[2] = im.vp_pos[2]*px + im.vp_pos[5]*py + im.vp_pos[11]
                for i in range(3):
                    ray.org[i] = v_pos[i]
                    ray.dir[i] = im.vp_dir[i]
                ray.tnear = 0.0
                ray.tfar = 1e37
                ray.geomID = rtcg.RTC_INVALID_GEOMETRY_ID
                ray.primID = rtcg.RTC_INVALID_GEOMETRY_ID
                ray.instID = rtcg.RTC_INVALID_GEOMETRY_ID
                ray.mask = -1
                ray.time = 0
                rtcs.rtcIntersect(scene.scene_i, ray)
                elemID = ray.primID / 12
                data[j] = ray.time
            self.aimage = data.reshape(self.image.nv[0], self.image.nv[1])
            free(v_pos)
        else:
            v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            # If we do not have a simple image plane, we have to cast all
            # our rays 
            for j in range(size):
                offset = j * 3
                for i in range(3): v_pos[i] = im.vp_pos[i + offset]
                for i in range(3): v_dir[i] = im.vp_dir[i + offset]
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
                rtcs.rtcIntersect(scene.scene_i, ray)
                data[j] = ray.time
            self.aimage = data.reshape(self.image.nv[0], self.image.nv[1])
            free(v_pos)
            free(v_dir)

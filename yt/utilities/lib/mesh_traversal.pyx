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


cdef void error_printer(const rtc.RTCError code, const char *_str):
    print "ERROR CAUGHT IN EMBREE"
    rtc.print_error(code)
    print "ERROR MESSAGE:", _str

cdef class EmbreeVolume:

    def __init__(self):
        rtc.rtcInit(NULL)
        rtc.rtcSetErrorFunction(error_printer)
        self.scene_i = rtcs.rtcNewScene(rtcs.RTC_SCENE_STATIC, rtcs.RTC_INTERSECT1)

    def __dealloc__(self):
        rtcs.rtcDeleteScene(self.scene_i)

cdef class MeshSampler(ImageSampler):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self, EmbreeVolume volume, int num_threads = 0):
        '''

        This function is supposed to cast the rays and return the
        image.

        '''

        rtcs.rtcCommit(volume.scene_i)
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, j, ni, nj, nn
        cdef np.int64_t offset
        cdef np.int64_t iter[4]
        cdef ImageContainer *im = self.image
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.float64_t rgba[6]
        cdef np.float64_t extrema[4]
        cdef np.float64_t max_t
        hit = 0
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
        cdef int vd_i = 0
        cdef int vd_step = 1
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
                offset = im.im_strides[0] * vi + im.im_strides[1] * vj
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
                vd_i += vd_step
                rtcs.rtcIntersect(volume.scene_i, ray)
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
                if v_dir[0] == v_dir[1] == v_dir[2] == 0.0:
                    continue

            free(v_dir)
            free(v_pos)
        return hit

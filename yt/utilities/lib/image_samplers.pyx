# distutils: include_dirs = LIB_DIR
# distutils: extra_compile_args = CPP14_FLAG OMP_ARGS
# distutils: extra_link_args = CPP14_FLAG OMP_ARGS
# distutils: libraries = STD_LIBS
# distutils: sources = FIXED_INTERP
# distutils: language = c++
"""
Image sampler definitions



"""


import numpy as np

cimport cython
cimport lenses
from field_interpolation_tables cimport (
    FieldInterpolationTable,
    FIT_eval_transfer,
    FIT_eval_transfer_with_light,
    FIT_initialize_table,
)
from libc.math cimport sqrt
from libc.stdlib cimport free, malloc

from yt.utilities.lib.fp_utils cimport fclip, i64clip, imin

from .fixed_interpolator cimport (
    eval_gradient,
    fast_interpolate,
    offset_fill,
    offset_interpolate,
    trilinear_interpolate,
    vertex_interp,
)
from .grid_traversal cimport sampler_function, walk_volume

from yt.funcs import mylog

from ._octree_raytracing cimport RayInfo, _OctreeRayTracing


cdef extern from "platform_dep.h":
    long int lrint(double x) nogil

DEF Nch = 4

from cython.parallel import parallel, prange, threadid

from cpython.exc cimport PyErr_CheckSignals
from vec3_ops cimport L2_norm, dot, fma, subtract


cdef struct VolumeRenderAccumulator:
    int n_fits
    int n_samples
    FieldInterpolationTable *fits
    int field_table_ids[6]
    np.float64_t star_coeff
    np.float64_t star_er
    np.float64_t star_sigma_num
    np.float64_t *light_dir
    np.float64_t *light_rgba
    int grey_opacity


cdef class ImageSampler:
    def __init__(self,
                  np.float64_t[:,:,:] vp_pos,
                  np.float64_t[:,:,:] vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  str volume_method,
                  *args,
                  **kwargs):
        cdef int i

        self.volume_method = volume_method
        camera_data = kwargs.pop("camera_data", None)
        if camera_data is not None:
            self.camera_data = camera_data

        zbuffer = kwargs.pop("zbuffer", None)
        if zbuffer is None:
            zbuffer = np.ones((image.shape[0], image.shape[1]), "float64")

        image_used = np.zeros((image.shape[0], image.shape[1]), "int64")
        mesh_lines = np.zeros((image.shape[0], image.shape[1]), "int64")

        self.lens_type = kwargs.pop("lens_type", None)
        if self.lens_type == "plane-parallel":
            self.extent_function = lenses.calculate_extent_plane_parallel
            self.vector_function = lenses.generate_vector_info_plane_parallel
        else:
            if not (vp_pos.shape[0] == vp_dir.shape[0] == image.shape[0]) or \
               not (vp_pos.shape[1] == vp_dir.shape[1] == image.shape[1]):
                msg = "Bad lens shape / direction for %s\n" % (self.lens_type)
                msg += "Shapes: (%s - %s - %s) and (%s - %s - %s)" % (
                    vp_pos.shape[0], vp_dir.shape[0], image.shape[0],
                    vp_pos.shape[1], vp_dir.shape[1], image.shape[1])
                raise RuntimeError(msg)

            if camera_data is not None and self.lens_type == 'perspective':
                self.extent_function = lenses.calculate_extent_perspective
            else:
                self.extent_function = lenses.calculate_extent_null
            self.\
                vector_function = lenses.generate_vector_info_null

        # These assignments are so we can track the objects and prevent their
        # de-allocation from reference counts.  Note that we do this to the
        # "atleast_3d" versions.  Also, note that we re-assign the input
        # arguments.
        self.vp_pos = vp_pos
        self.vp_dir = vp_dir
        self.image = self.aimage = image
        self.acenter = center
        self.center = <np.float64_t *> center.data
        self.ax_vec = x_vec
        self.x_vec = <np.float64_t *> x_vec.data
        self.ay_vec = y_vec
        self.y_vec = <np.float64_t *> y_vec.data
        self.zbuffer = zbuffer
        self.azbuffer = np.asarray(zbuffer)
        self.image_used = image_used
        self.aimage_used = np.asarray(image_used)
        self.mesh_lines = mesh_lines
        self.amesh_lines = np.asarray(mesh_lines)
        self.nv[0] = image.shape[0]
        self.nv[1] = image.shape[1]
        for i in range(4): self.bounds[i] = bounds[i]
        self.pdx = (bounds[1] - bounds[0])/self.nv[0]
        self.pdy = (bounds[3] - bounds[2])/self.nv[1]
        for i in range(3):
            self.width[i] = width[i]

    def __call__(self, PartitionedGrid pg, **kwa):
        if self.volume_method == 'KDTree':
            return self.cast_through_kdtree(pg, **kwa)
        elif self.volume_method == 'Octree':
            return self.cast_through_octree(pg, **kwa)
        else:
            raise NotImplementedError(
                'Volume rendering has not been implemented for method: "%s"' %
                self.volume_method
            )

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def cast_through_kdtree(self, PartitionedGrid pg, int num_threads = 0):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, j
        cdef VolumeContainer *vc = pg.container
        self.setup(pg)
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.float64_t max_t
        hit = 0
        cdef np.int64_t nx, ny, size
        cdef np.int64_t iter[4]
        self.extent_function(self, vc, iter)
        iter[0] = i64clip(iter[0]-1, 0, self.nv[0])
        iter[1] = i64clip(iter[1]+1, 0, self.nv[0])
        iter[2] = i64clip(iter[2]-1, 0, self.nv[1])
        iter[3] = i64clip(iter[3]+1, 0, self.nv[1])
        nx = (iter[1] - iter[0])
        ny = (iter[3] - iter[2])
        size = nx * ny
        cdef ImageAccumulator *idata
        cdef np.float64_t width[3]
        cdef int chunksize = 100
        for i in range(3):
            width[i] = self.width[i]
        with nogil, parallel(num_threads = num_threads):
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            idata.supp_data = self.supp_data
            v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            for j in prange(size, schedule="static", chunksize=chunksize):
                vj = j % ny
                vi = (j - vj) / ny + iter[0]
                vj = vj + iter[2]
                # Dynamically calculate the position
                self.vector_function(self, vi, vj, width, v_dir, v_pos)
                for i in range(Nch):
                    idata.rgba[i] = self.image[vi, vj, i]
                max_t = fclip(self.zbuffer[vi, vj], 0.0, 1.0)
                walk_volume(vc, v_pos, v_dir, self.sample,
                            (<void *> idata), NULL, max_t)
                if (j % (10*chunksize)) == 0:
                    with gil:
                        PyErr_CheckSignals()
                for i in range(Nch):
                    self.image[vi, vj, i] = idata.rgba[i]
            idata.supp_data = NULL
            free(idata)
            free(v_pos)
            free(v_dir)
        return hit

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def cast_through_octree(self, PartitionedGrid pg, _OctreeRayTracing oct, int num_threads = 0):
        cdef RayInfo[int]* ri
        self.setup(pg)

        cdef sampler_function* sampler = <sampler_function*> self.sample

        cdef np.int64_t nx, ny, size
        cdef VolumeContainer *vc
        cdef ImageAccumulator *idata

        nx = self.nv[0]
        ny = self.nv[1]
        size = nx * ny

        cdef int i, j, k, vi, vj, icell
        cdef int[3] index = [0, 0, 0]
        cdef int chunksize = 100

        cdef int n_fields = pg.container.n_fields

        cdef np.float64_t vp_dir_len

        mylog.debug('Integrating rays')
        with nogil, parallel(num_threads=num_threads):
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            idata.supp_data = self.supp_data

            ri = new RayInfo[int]()

            vc = <VolumeContainer*> malloc(sizeof(VolumeContainer))
            vc.n_fields = 1
            vc.data = <np.float64_t**> malloc(sizeof(np.float64_t*))
            vc.mask = <np.uint8_t*> malloc(8*sizeof(np.uint8_t))
            # The actual dimensions are 2x2x2, but the sampler
            # assumes vertex-centred data for a 1x1x1 lattice (i.e.
            # 2^3 vertices)
            vc.dims[0] = 1
            vc.dims[1] = 1
            vc.dims[2] = 1
            for j in prange(size, schedule='static', chunksize=chunksize):
                vj = j % ny
                vi = (j - vj) / ny

                vp_dir_len = sqrt(
                    self.vp_dir[vi, vj, 0]**2 +
                    self.vp_dir[vi, vj, 1]**2 +
                    self.vp_dir[vi, vj, 2]**2)

                # Cast ray
                oct.oct.cast_ray(&self.vp_pos[vi, vj, 0], &self.vp_dir[vi, vj, 0],
                                 ri.keys, ri.t)
                # Contains the ordered indices of the cells hit by the ray
                # and the entry/exit t values
                if ri.keys.size() == 0:
                    continue

                for i in range(Nch):
                    idata.rgba[i] = 0
                for i in range(8):
                    vc.mask[i] = 1

                # Iterate over cells
                for i in range(ri.keys.size()):
                    icell = ri.keys[i]
                    for k in range(n_fields):
                        vc.data[k] = &pg.container.data[k][14*icell]

                    # Fill the volume container with the current boundaries
                    for k in range(3):
                        vc.left_edge[k] = pg.container.data[0][14*icell+8+k]
                        vc.right_edge[k] = pg.container.data[0][14*icell+11+k]
                        vc.dds[k] = (vc.right_edge[k] - vc.left_edge[k])
                        vc.idds[k] = 1/vc.dds[k]
                    # Now call the sampler
                    sampler(
                        vc,
                        &self.vp_pos[vi, vj, 0],
                        &self.vp_dir[vi, vj, 0],
                        ri.t[2*i  ]/vp_dir_len,
                        ri.t[2*i+1]/vp_dir_len,
                        index,
                        <void *> idata
                        )
                for i in range(Nch):
                    self.image[vi, vj, i] = idata.rgba[i]

                # Empty keys and t
                ri.keys.clear()
                ri.t.clear()

            del ri
            free(vc.data)
            free(vc.mask)
            free(vc)
            idata.supp_data = NULL
            free(idata)

        mylog.debug('Done integration')


    cdef void setup(self, PartitionedGrid pg):
        return

    @staticmethod
    cdef void sample(
                 VolumeContainer *vc,
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
        return

    def ensure_code_unit_params(self, params):
        for param_name in ['center', 'vp_pos', 'vp_dir', 'width']:
            param = params[param_name]
            if hasattr(param, 'in_units'):
                params[param_name] = param.in_units('code_length')
        bounds = params['bounds']
        if hasattr(bounds[0], 'units'):
            params['bounds'] = tuple(b.in_units('code_length').d for b in bounds)

        return params

cdef class ProjectionSampler(ImageSampler):

    @staticmethod
    cdef void sample(
                 VolumeContainer *vc,
                 np.float64_t v_pos[3],
                 np.float64_t v_dir[3],
                 np.float64_t enter_t,
                 np.float64_t exit_t,
                 int index[3],
                 void *data) nogil:
        cdef ImageAccumulator *im = <ImageAccumulator *> data
        cdef int i
        cdef np.float64_t dl = (exit_t - enter_t)
        cdef int di = (index[0]*vc.dims[1]+index[1])*vc.dims[2]+index[2]
        for i in range(imin(4, vc.n_fields)):
            im.rgba[i] += vc.data[i][di] * dl


cdef class InterpolatedProjectionSampler(ImageSampler):
    def __cinit__(self,
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  str volume_method,
                  n_samples = 10,
                  **kwargs
        ):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, volume_method, **kwargs)
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.n_samples = n_samples
        self.supp_data = <void *> self.vra

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @staticmethod
    cdef void sample(
                     VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     np.float64_t enter_t,
                     np.float64_t exit_t,
                     int index[3],
                     void *data) nogil:
        cdef ImageAccumulator *im = <ImageAccumulator *> data
        cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
                im.supp_data
        # we assume this has vertex-centered data.
        cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                        + index[1] * (vc.dims[2] + 1) + index[2]
        cdef np.float64_t dp[3]
        cdef np.float64_t ds[3]
        cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
        cdef np.float64_t dvs[6]
        for i in range(3):
            dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
            dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
            dp[i] *= vc.idds[i]
            ds[i] = v_dir[i] * vc.idds[i] * dt
        for i in range(vri.n_samples):
            for j in range(vc.n_fields):
                dvs[j] = offset_interpolate(vc.dims, dp,
                        vc.data[j] + offset)
            for j in range(imin(3, vc.n_fields)):
                im.rgba[j] += dvs[j] * dt
            for j in range(3):
                dp[j] += ds[j]


cdef class VolumeRenderSampler(ImageSampler):
    def __cinit__(self,
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  str volume_method,
                  tf_obj,
                  n_samples = 10,
                  **kwargs
        ):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, volume_method, **kwargs)
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] temp
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.fits = <FieldInterpolationTable *> \
            malloc(sizeof(FieldInterpolationTable) * 6)
        self.vra.n_fits = tf_obj.n_field_tables
        assert(self.vra.n_fits <= 6)
        self.vra.grey_opacity = getattr(tf_obj, "grey_opacity", 0)
        self.vra.n_samples = n_samples
        self.my_field_tables = []
        for i in range(self.vra.n_fits):
            temp = tf_obj.tables[i].y
            FIT_initialize_table(&self.vra.fits[i],
                      temp.shape[0],
                      <np.float64_t *> temp.data,
                      tf_obj.tables[i].x_bounds[0],
                      tf_obj.tables[i].x_bounds[1],
                      tf_obj.field_ids[i], tf_obj.weight_field_ids[i],
                      tf_obj.weight_table_ids[i])
            self.my_field_tables.append((tf_obj.tables[i],
                                         tf_obj.tables[i].y))
        for i in range(6):
            self.vra.field_table_ids[i] = tf_obj.field_table_ids[i]
        self.supp_data = <void *> self.vra

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @staticmethod
    cdef void sample(
                     VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     np.float64_t enter_t,
                     np.float64_t exit_t,
                     int index[3],
                     void *data) nogil:
        cdef ImageAccumulator *im = <ImageAccumulator *> data
        cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
                im.supp_data
        # we assume this has vertex-centered data.
        cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                        + index[1] * (vc.dims[2] + 1) + index[2]
        cdef int cell_offset = index[0] * (vc.dims[1]) * (vc.dims[2]) \
                        + index[1] * (vc.dims[2]) + index[2]
        if vc.mask[cell_offset] != 1:
            return
        cdef np.float64_t dp[3]
        cdef np.float64_t ds[3]
        cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
        cdef np.float64_t dvs[6]
        for i in range(3):
            dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
            dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
            dp[i] *= vc.idds[i]
            ds[i] = v_dir[i] * vc.idds[i] * dt
        for i in range(vri.n_samples):
            for j in range(vc.n_fields):
                dvs[j] = offset_interpolate(vc.dims, dp,
                        vc.data[j] + offset)
            FIT_eval_transfer(dt, dvs, im.rgba, vri.n_fits,
                    vri.fits, vri.field_table_ids, vri.grey_opacity)
            for j in range(3):
                dp[j] += ds[j]

    def __dealloc__(self):
        for i in range(self.vra.n_fits):
            free(self.vra.fits[i].d0)
            free(self.vra.fits[i].dy)
        free(self.vra.fits)
        free(self.vra)

cdef class LightSourceRenderSampler(ImageSampler):
    def __cinit__(self,
                  np.ndarray vp_pos,
                  np.ndarray vp_dir,
                  np.ndarray[np.float64_t, ndim=1] center,
                  bounds,
                  np.ndarray[np.float64_t, ndim=3] image,
                  np.ndarray[np.float64_t, ndim=1] x_vec,
                  np.ndarray[np.float64_t, ndim=1] y_vec,
                  np.ndarray[np.float64_t, ndim=1] width,
                  str volume_method,
                  tf_obj,
                  n_samples = 10,
                  light_dir=[1.,1.,1.],
                  light_rgba=[1.,1.,1.,1.],
                  **kwargs):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, volume_method, **kwargs)
        cdef int i
        cdef np.ndarray[np.float64_t, ndim=1] temp
        # Now we handle tf_obj
        self.vra = <VolumeRenderAccumulator *> \
            malloc(sizeof(VolumeRenderAccumulator))
        self.vra.fits = <FieldInterpolationTable *> \
            malloc(sizeof(FieldInterpolationTable) * 6)
        self.vra.n_fits = tf_obj.n_field_tables
        assert(self.vra.n_fits <= 6)
        self.vra.grey_opacity = getattr(tf_obj, "grey_opacity", 0)
        self.vra.n_samples = n_samples
        self.vra.light_dir = <np.float64_t *> malloc(sizeof(np.float64_t) * 3)
        self.vra.light_rgba = <np.float64_t *> malloc(sizeof(np.float64_t) * 4)
        light_dir /= np.sqrt(light_dir[0] * light_dir[0] +
                             light_dir[1] * light_dir[1] +
                             light_dir[2] * light_dir[2])
        for i in range(3):
            self.vra.light_dir[i] = light_dir[i]
        for i in range(4):
            self.vra.light_rgba[i] = light_rgba[i]
        self.my_field_tables = []
        for i in range(self.vra.n_fits):
            temp = tf_obj.tables[i].y
            FIT_initialize_table(&self.vra.fits[i],
                      temp.shape[0],
                      <np.float64_t *> temp.data,
                      tf_obj.tables[i].x_bounds[0],
                      tf_obj.tables[i].x_bounds[1],
                      tf_obj.field_ids[i], tf_obj.weight_field_ids[i],
                      tf_obj.weight_table_ids[i])
            self.my_field_tables.append((tf_obj.tables[i],
                                         tf_obj.tables[i].y))
        for i in range(6):
            self.vra.field_table_ids[i] = tf_obj.field_table_ids[i]
        self.supp_data = <void *> self.vra

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @staticmethod
    cdef void sample(
                     VolumeContainer *vc,
                     np.float64_t v_pos[3],
                     np.float64_t v_dir[3],
                     np.float64_t enter_t,
                     np.float64_t exit_t,
                     int index[3],
                     void *data) nogil:
        cdef ImageAccumulator *im = <ImageAccumulator *> data
        cdef VolumeRenderAccumulator *vri = <VolumeRenderAccumulator *> \
                im.supp_data
        # we assume this has vertex-centered data.
        cdef int offset = index[0] * (vc.dims[1] + 1) * (vc.dims[2] + 1) \
                        + index[1] * (vc.dims[2] + 1) + index[2]
        cdef np.float64_t dp[3]
        cdef np.float64_t ds[3]
        cdef np.float64_t dt = (exit_t - enter_t) / vri.n_samples
        cdef np.float64_t dvs[6]
        cdef np.float64_t *grad
        grad = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
        for i in range(3):
            dp[i] = (enter_t + 0.5 * dt) * v_dir[i] + v_pos[i]
            dp[i] -= index[i] * vc.dds[i] + vc.left_edge[i]
            dp[i] *= vc.idds[i]
            ds[i] = v_dir[i] * vc.idds[i] * dt
        for i in range(vri.n_samples):
            for j in range(vc.n_fields):
                dvs[j] = offset_interpolate(vc.dims, dp,
                        vc.data[j] + offset)
            eval_gradient(vc.dims, dp, vc.data[0] + offset, grad)
            FIT_eval_transfer_with_light(dt, dvs, grad,
                    vri.light_dir, vri.light_rgba,
                    im.rgba, vri.n_fits,
                    vri.fits, vri.field_table_ids, vri.grey_opacity)
            for j in range(3):
                dp[j] += ds[j]
        free(grad)


    def __dealloc__(self):
        for i in range(self.vra.n_fits):
            free(self.vra.fits[i].d0)
            free(self.vra.fits[i].dy)
        free(self.vra.light_dir)
        free(self.vra.light_rgba)
        free(self.vra.fits)
        free(self.vra)

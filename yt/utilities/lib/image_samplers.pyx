"""
Image sampler definitions



"""


import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, calloc, free, abs
from libc.math cimport exp, floor, log2, \
    fabs, atan, atan2, asin, cos, sin, sqrt, acos, M_PI
from libcpp.vector cimport vector
from yt.utilities.lib.fp_utils cimport imax, fmax, imin, fmin, iclip, fclip, i64clip
from yt.utilities.lib.octree_raytracing cimport \
    ray_step, DomainDecomposition_find_domain, OctreesDescriptor, OctreeDescriptor
from yt.geometry.oct_container cimport SparseOctreeContainer
from yt.frontends.ramses.io_utils cimport hilbert3d_single

from field_interpolation_tables cimport \
    FieldInterpolationTable, FIT_initialize_table, FIT_eval_transfer,\
    FIT_eval_transfer_with_light
cimport lenses
from .grid_traversal cimport walk_volume, sampler_function
from .fixed_interpolator cimport \
    offset_interpolate, \
    fast_interpolate, \
    trilinear_interpolate, \
    eval_gradient, \
    offset_fill, \
    vertex_interp

cdef extern from "platform_dep.h":
    long int lrint(double x) nogil

DEF Nch = 4

from cython.parallel import prange, parallel, threadid
from vec3_ops cimport dot, subtract, L2_norm, fma

from cpython.exc cimport PyErr_CheckSignals

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

# Dummy function to call static methods from Cython
cdef void walk_cell(VolumeContainer *vc, np.float64_t v_pos[3], np.float64_t v_dir[3],
                    np.float64_t enter_t, np.float64_t exit_t,
                    int index[3],
                    sampler_function *sample, void *data) nogil:
    sample(vc, v_pos, v_dir, enter_t, exit_t, index, data)

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
                  *args, **kwargs):
        cdef int i

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def __call__(self, PartitionedGrid pg, int num_threads = 0):
        # This routine will iterate over all of the vectors and cast each in
        # turn.  Might benefit from a more sophisticated intersection check,
        # like http://courses.csusm.edu/cs697exz/ray_box.htm
        cdef int vi, vj, hit, i, j
        cdef np.int64_t iter[4]
        cdef VolumeContainer *vc = pg.container
        self.setup(pg)
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.float64_t max_t
        hit = 0
        cdef np.int64_t nx, ny, size
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
    def call_octrees(self, data_source, OctreesDescriptor octrees_descriptor, int num_threads = 0):
        cdef int vi, vj, hit, i, j
        cdef np.int64_t iter[4]
        cdef VolumeContainer *vc = NULL
        # self.setup(pg) # TODO
        cdef np.float64_t *v_pos
        cdef np.float64_t *v_dir
        cdef np.float64_t max_t
        hit = 0
        cdef np.int64_t nx, ny, size
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

        # Domain decomposition information
        cdef np.float64_t bscale
        cdef int bit_length
        cdef int ncpu
        cdef np.uint64_t[:] keys
        cdef SparseOctreeContainer octree
        cdef np.float64_t[3] DLE, DRE
        # Get the first one for DLE and DRE, this will be overriden below
        self.octrees_descriptor = octrees_descriptor
        octree = <SparseOctreeContainer> self.octrees_descriptor.octrees[1]

        DLE = octree.DLE
        DRE = octree.DRE
        ds = data_source.ds
        bscale = ds.hilbert['bscale']
        bit_length = ds.hilbert['bit_length']
        ncpu = ds.parameters['ncpu']
        keys = ds.hilbert['keys'].astype(np.uint64)

        # Containers
        cdef vector[np.uint64_t] octList, countPerDomain, domainList
        cdef vector[np.uint8_t] cellList
        cdef vector[np.float64_t] tList

        for i in range(3):
            width[i] = self.width[i]
        with nogil, parallel(num_threads = num_threads):
            idata = <ImageAccumulator *> malloc(sizeof(ImageAccumulator))
            idata.supp_data = self.supp_data
            v_pos = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            v_dir = <np.float64_t *> malloc(3 * sizeof(np.float64_t))
            vc = <VolumeContainer *> malloc(sizeof(VolumeContainer))
            vc.n_fields = 1 # FIXME
            vc.mask = <np.uint8_t *> malloc(sizeof(np.uint8_t)*8)
            vc.data = <np.float64_t **> malloc(sizeof(np.float64_t*))      # only 1 field is supported
            vc.data[0] = <np.float64_t *> malloc(sizeof(np.float64_t) * 27)# 3x3x3
            vc.dims[0] = 2
            vc.dims[1] = 2
            vc.dims[2] = 2
            for j in prange(size, schedule="static", chunksize=chunksize):
                vj = j % ny
                vi = (j - vj) / ny + iter[0]
                vj = vj + iter[2]
                # Dynamically calculate the position
                self.vector_function(self, vi, vj, width, v_dir, v_pos)
                for i in range(Nch):
                    idata.rgba[i] = self.image[vi, vj, i]
                max_t = fclip(self.zbuffer[vi, vj], 0.0, 1.0)

                self.octree_cast_single_ray(
                    vc, octList, countPerDomain, domainList, cellList, tList,
                    v_pos, v_dir, DLE, DRE, bscale, bit_length, ncpu, keys,
                    (<void *> idata))


                if (j % (10*chunksize)) == 0:
                    with gil:
                        PyErr_CheckSignals()
                for i in range(Nch):
                    self.image[vi, vj, i] = idata.rgba[i]
            idata.supp_data = NULL
            free(idata)
            free(v_pos)
            free(v_dir)
            free(vc.mask)
            free(vc.data)
            free(vc)

        del self.octrees
        del self.octree_descriptors

        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void octree_cast_single_ray(
            self,
            VolumeContainer* vc,
            vector[np.uint64_t] octList,
            vector[np.uint64_t] countPerDomain,
            vector[np.uint64_t] domainList,
            vector[np.uint8_t] cellList,
            vector[np.float64_t] tList,
            np.float64_t* v_pos,
            np.float64_t* v_dir,
            np.float64_t* DLE,
            np.float64_t* DRE,
            np.float64_t bscale,
            int bit_length,
            int ncpu,
            np.uint64_t[:] keys,
            void *idata
            ) nogil:

        cdef int nAdded, count, nextDom
        cdef np.float64_t tmin, tin
        cdef np.float64_t[3] pos
        cdef int i, j, k, l, iglobal

        # Clear container
        octList.clear()
        cellList.clear()
        tList.clear()
        countPerDomain.clear()
        domainList.clear()

        # Find entrance hit node
        tin = 0
        for i in range(3):
            if v_dir[i] < 0:
                tin = max(tin, DRE[i] - v_pos[i]) / v_dir[i]
            else:
                tin = max(tin, DLE[i] - v_pos[i]) / v_dir[i]
            pos[i] = v_pos[i] + tin * v_dir[i]

        with gil:
            # For some reason, cannot do that w/o the gil...
            nextDom = DomainDecomposition_find_domain(bscale, bit_length, ncpu, keys, pos)

        tmin = 0
        count = 0
        nAdded = 0
        while nextDom > 0:
            # Add domain to list of domains
            domainList.push_back(nextDom)

            # Call ray traversal on domain
            with gil:
                octree = <SparseOctreeContainer> self.octrees[nextDom]
            nextDom = ray_step(octree, v_pos, v_dir, octList, cellList, tList, tmin, nextDom)

            # Update number of cells crossed
            nAdded = tList.size() - count
            count = tList.size()

            # Store this number
            countPerDomain.push_back(nAdded)

            # Next starting point
            if tList.size() > 0:
                tmin = tList.back()

        cdef int ioct, prev_ioct, idom
        cdef int ind[3]
        cdef OctreeDescriptor *od
        cdef np.float64_t[:, :] data
        iglobal = 0
        for i in range(domainList.size()):
            idom = domainList[i]
            # Get data
            od = self.octrees_descriptor.get_descriptor(idom)

            prev_ioct = -1
            # Loop over cells
            for j in range(countPerDomain[i]):
                ioct = octList[iglobal]

                if ioct != prev_ioct and ioct > -1:
                    prev_ioct = ioct

                    for k in range(3):
                        vc.left_edge[k] = od.oct_LE[ioct*64+k]
                        vc.dds[k] = od.oct_ddd[ioct*64+k]
                        vc.idds[k] = 1/vc.dds[k]

                    for k in range(od.nfields):
                        vc.data[k] = &od.data[ioct*64+k]

                ind[0] = (cellList[iglobal] & 0b100) >> 2
                ind[1] = (cellList[iglobal] & 0b010) >> 1
                ind[2] = (cellList[iglobal] & 0b001) >> 0

                walk_cell(vc, v_pos, v_dir, tList[iglobal], tList[iglobal+1], ind,
                          self.sample, idata)

                iglobal += 1

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
                  n_samples = 10, **kwargs):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, **kwargs)
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
                  tf_obj, n_samples = 10,
                  **kwargs):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, **kwargs)
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
                  tf_obj, n_samples = 10,
                  light_dir=[1.,1.,1.],
                  light_rgba=[1.,1.,1.,1.],
                  **kwargs):
        ImageSampler.__init__(self, vp_pos, vp_dir, center, bounds, image,
                               x_vec, y_vec, width, **kwargs)
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

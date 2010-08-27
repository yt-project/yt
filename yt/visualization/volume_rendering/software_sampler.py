"""
Import the components of the volume rendering extension

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py
import numpy as na

from yt.funcs import *

from yt.data_objects.data_containers import data_object_registry
from yt.utilities.amr_utils import TransferFunctionProxy, VectorPlane
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.visualization.volume_rendering.grid_partitioner import \
    HomogenizedBrickCollection

# We're going to register this class, but it does not directly inherit from
# AMRData.
class VolumeRendering(ParallelAnalysisInterface):
    bricks = None
    def __init__(self, normal_vector, width, center,
                 resolution, transfer_function,
                 fields = None, whole_box = False,
                 sub_samples = 5, north_vector = None,
                 pf = None):
        # Now we replicate some of the 'cutting plane' logic
        if not iterable(resolution):
            resolution = (resolution, resolution)
        self.resolution = resolution
        self.sub_samples = sub_samples
        if not iterable(width):
            width = (width, width, width) # front/back, left/right, top/bottom
        self.width = width
        self.center = center
        if fields is None: fields = ["Density"]
        self.fields = fields
        self.transfer_function = transfer_function

        # Now we set up our  various vectors
        normal_vector /= na.sqrt( na.dot(normal_vector, normal_vector))
        if north_vector is None:
            vecs = na.identity(3)
            t = na.cross(normal_vector, vecs).sum(axis=1)
            ax = t.argmax()
            north_vector = na.cross(vecs[ax,:], normal_vector).ravel()
        north_vector /= na.sqrt(na.dot(north_vector, north_vector))
        east_vector = -na.cross(north_vector, normal_vector).ravel()
        east_vector /= na.sqrt(na.dot(east_vector, east_vector))
        self.unit_vectors = [north_vector, east_vector, normal_vector]
        self.box_vectors = na.array([self.unit_vectors[0]*self.width[0],
                                     self.unit_vectors[1]*self.width[1],
                                     self.unit_vectors[2]*self.width[2]])

        self.origin = center - 0.5*width[0]*self.unit_vectors[0] \
                             - 0.5*width[1]*self.unit_vectors[1] \
                             - 0.5*width[2]*self.unit_vectors[2]
        self.back_center = center - 0.5*width[0]*self.unit_vectors[2]
        self.front_center = center + 0.5*width[0]*self.unit_vectors[2]

        self._initialize_source()
        self._construct_vector_array()

    def _initialize_source(self):
        check, source, rf = self._partition_hierarchy_2d_inclined(
                self.unit_vectors, self.origin, self.width, self.box_vectors)
        if check:
            self._base_source = self.pf.h.inclined_box(
                self.origin, self.box_vectors)
        else:
            # To avoid doubling-up
            self._base_source = source
        self.source = source
        self.res_fac = rf
        # Note that if we want to do this in parallel, with 3D domain decomp
        # for the grid/bricks, we can supply self._base_source here.  But,
        # _distributed can't be overridden in that case.
        self._brick_collection = HomogenizedBrickCollection(self.source)

    def ray_cast(self, finalize=True):
        if self.bricks is None: self.partition_grids()
        # Now we order our bricks
        total_cells, LE, RE = 0, [], []
        for b in self.bricks:
            LE.append(b.LeftEdge)
            RE.append(b.RightEdge)
            total_cells += na.prod(b.my_data[0].shape)
        LE = na.array(LE) - self.back_center
        RE = na.array(RE) - self.back_center
        LE = na.sum(LE * self.unit_vectors[2], axis=1)
        RE = na.sum(RE * self.unit_vectors[2], axis=1)
        dist = na.minimum(LE, RE)
        ind = na.argsort(dist)
        pbar = get_pbar("Ray casting ", total_cells)
        total_cells = 0
        tfp = TransferFunctionProxy(self.transfer_function)
        tfp.ns = self.sub_samples
        for i, b in enumerate(self.bricks[ind]):
            pos = b.cast_plane(tfp, self.vector_plane)
            total_cells += na.prod(b.my_data[0].shape)
            pbar.update(total_cells)
        pbar.finish()
        if finalize: self._finalize()

    def _finalize(self):
        #im = self._mpi_catdict(dict(image=self.image)).pop('image')
        im, f = self._mpi_catrgb((self.image, self.resolution))
        self.image = im

    def dump_image(self, prefix):
        fn = "%s.h5" % (self._get_filename(prefix))
        mylog.info("Saving to %s", fn)
        f = h5py.File(fn, "w")
        f.create_dataset("/image", data=self.image)

    def load_bricks(self, fn):
        self.bricks = import_partitioned_grids(fn)

    def save_bricks(self, fn):
        # This will need to be modified for parallel
        export_partitioned_grids(self.bricks, fn)

    def save_image(self, prefix = None, norm = 1.0):
        if norm is not None:
            mi, ma = self.image.min(), norm*self.image.max()
            print "Normalizing with ", mi, ma
            image = (na.clip(self.image, mi, ma) - mi)/(ma - mi)
        else:
            image = self.image
        if prefix is None: prefix = "%s_volume_rendering" % (self.pf)
        plot_rgb(image, prefix)

    def partition_grids(self):
        log_field = []
        for field in self.fields:
            log_field.append(field in self.pf.field_info and 
                             self.pf.field_info[field].take_log)
        self._brick_collection._partition_local_grids(self.fields, log_field)
        # UNCOMMENT FOR PARALLELISM
        #self._brick_collection._collect_bricks(self.source)
        self.bricks = self._brick_collection.bricks

    def _construct_vector_array(self):
        rx = self.resolution[0] * self.res_fac[0]
        ry = self.resolution[1] * self.res_fac[1]
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        self.image = na.zeros((rx,ry,3), dtype='float64', order='C')
        # We might have a different width and back_center
        bl = self.source.box_lengths
        px = na.linspace(-bl[0]/2.0, bl[0]/2.0, rx)[:,None]
        py = na.linspace(-bl[1]/2.0, bl[1]/2.0, ry)[None,:]
        inv_mat = self.source._inv_mat
        bc = self.source.origin + 0.5*self.source.box_vectors[0] \
                                + 0.5*self.source.box_vectors[1]
        vectors = na.zeros((rx, ry, 3),
                            dtype='float64', order='C')
        vectors[:,:,0] = inv_mat[0,0]*px + inv_mat[0,1]*py + bc[0]
        vectors[:,:,1] = inv_mat[1,0]*px + inv_mat[1,1]*py + bc[1]
        vectors[:,:,2] = inv_mat[2,0]*px + inv_mat[2,1]*py + bc[2]
        bounds = (px.min(), px.max(), py.min(), py.max())
        self.vector_plane = VectorPlane(vectors, self.box_vectors[2],
                                    bc, bounds, self.image,
                                    self.source._x_vec, self.source._y_vec)
        self.vp_bounds = bounds
        self.vectors = vectors

data_object_registry["volume_rendering"] = VolumeRendering

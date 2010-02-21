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

import numpy as na
from yt.extensions.volume_rendering import *
from yt.funcs import *
from yt.lagos import data_object_registry, ParallelAnalysisInterface

def direct_ray_cast(pf, L, center, W, Nvec, tf, 
                    partitioned_grids = None, field = 'Density',
                    log_field = True, whole_box=False, region=None,
                    nsamples = 5):
    center = na.array(center, dtype='float64')

    # This just helps us keep track of stuff, and it's cheap
    cp = pf.h.cutting(L, center)
    back_center = center - cp._norm_vec * na.sqrt(3) * W
    front_center = center + cp._norm_vec * na.sqrt(3) *  W
    if region is None:
        if whole_box:
            cylinder = pf.h.region([0.5]*3,[0.0]*3,[1.0]*3)
        else:
            cylinder = pf.h.disk(center, L, na.sqrt(3)*W, 2*W*na.sqrt(3))
    else:
        if whole_box:
            print 'Warning, using region that you specified, not the whole box'
        cylinder = pf.h.region(region[0],region[1],region[2])

    if partitioned_grids == None:
        partitioned_grids = partition_all_grids(cylinder._grids,
                                    eval_func = lambda g: na.any(cylinder._get_point_indices(g)),
                                                field = field, log_field = log_field)
    #partitioned_grids = partition_all_grids(pf.h.grids)

    LE = (na.array([grid.LeftEdge for grid in partitioned_grids]) - back_center) * cp._norm_vec
    RE = (na.array([grid.RightEdge for grid in partitioned_grids]) - back_center) * cp._norm_vec
    DL = na.sum(LE, axis=1); del LE
    DR = na.sum(RE, axis=1); del RE
    dist = na.minimum(DL, DR)
    ind = na.argsort(dist)
    
    image = na.zeros((Nvec,Nvec,4), dtype='float64', order='F')
    image[:,:,3] = 0.0

    # Now we need to generate regular x,y,z values in regular space for our vector
    # starting places.
    px, py = na.mgrid[-W:W:Nvec*1j,-W:W:Nvec*1j]
    xv = cp._inv_mat[0,0]*px + cp._inv_mat[0,1]*py + back_center[0]
    yv = cp._inv_mat[1,0]*px + cp._inv_mat[1,1]*py + back_center[1]
    zv = cp._inv_mat[2,0]*px + cp._inv_mat[2,1]*py + back_center[2]
    vectors = na.array([xv, yv, zv], dtype='float64').transpose()
    vectors = vectors.copy('F')
    xp0, xp1 = px.min(), px.max()
    yp0, yp1 = py.min(), py.max()

    ng = partitioned_grids.size
    norm_vec = cp._norm_vec
    norm_vec = cp._norm_vec * (2.0*W*na.sqrt(3))
    hit = 0
    tnow = time.time()

    vp = VectorPlane(vectors, norm_vec, back_center,
                     (xp0, xp1, yp0, yp1), image, cp._x_vec, cp._y_vec)

    tf.light_dir = cp._norm_vec + 0.5 * cp._x_vec + 0.5 * cp._y_vec
    cx, cy, cz = 0.3, -0.3, 0.3
    tf.light_dir = (cp._inv_mat[0,0]*cx + cp._inv_mat[0,1]*cy + cz,
                    cp._inv_mat[1,0]*cx + cp._inv_mat[1,1]*cy + cz,
                    cp._inv_mat[2,0]*cx + cp._inv_mat[2,1]*cy + cz)
    
    tfp = TransferFunctionProxy(tf)
    tfp.ns = nsamples


    total_cells = sum(na.prod(g.my_data.shape) for g in partitioned_grids)
    pbar = get_pbar("Ray casting ", total_cells)
    total_cells = 0
    for i,g in enumerate(partitioned_grids[ind]):
        pbar.update(total_cells)
        pos = g.cast_plane(tfp, vp)
        total_cells += na.prod(g.my_data.shape)
    pbar.finish()

    return partitioned_grids[ind], image, vectors, norm_vec, pos

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
        self.source = source
        self.res_fac = rf

    def ray_cast(self):
        if self.bricks is None: self.partition_grids()
        # Now we order our bricks
        total_cells, LE, RE = 0, [], []
        for b in self.bricks:
            LE.append(b.LeftEdge)
            RE.append(b.RightEdge)
            total_cells += na.prod(b.my_data.shape)
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
            total_cells += na.prod(b.my_data.shape)
            pbar.update(total_cells)
        pbar.finish()
        self._finalize()

    def _finalize(self):
        #im = self._mpi_catdict(dict(image=self.image)).pop('image')
        im, f = self._mpi_catrgb((self.image, self.resolution))
        self.image = im

    def load_bricks(self, fn):
        self.bricks = import_partitioned_grids(fn)

    def save_bricks(self, fn):
        # This will need to be modified for parallel
        export_partitioned_grids(self.bricks, fn)

    def partition_grids(self):
        log_field = (self.fields[0] in self.pf.field_info and 
                     self.pf.field_info[self.fields[0]].take_log)
        self.bricks = partition_all_grids(self.source._grids,
                            field = self.fields[0],
                            log_field = log_field)

    def _construct_vector_array(self):
        rx = self.resolution[0] * self.res_fac[0]
        ry = self.resolution[1] * self.res_fac[1]
        # We should move away from pre-generation of vectors like this and into
        # the usage of on-the-fly generation in the VolumeIntegrator module
        self.image = na.zeros((ry,rx,4), dtype='float64', order='F')
        # We might have a different width and back_center
        bl = self.source.box_lengths
        px = na.linspace(-bl[0]/2.0, bl[0]/2.0, rx)
        py = na.linspace(-bl[1]/2.0, bl[1]/2.0, ry)
        inv_mat = self.source._inv_mat
        bc = self.source.origin + 0.5*self.source.box_vectors[0] \
                                + 0.5*self.source.box_vectors[1]
        vectors = na.zeros((ry, rx, 3),
                            dtype='float64', order='F')
        vectors[:,:,0] = inv_mat[0,0]*px[None,:] + inv_mat[0,1]*py[:,None] + bc[0]
        vectors[:,:,1] = inv_mat[1,0]*px[None,:] + inv_mat[1,1]*py[:,None] + bc[1]
        vectors[:,:,2] = inv_mat[2,0]*px[None,:] + inv_mat[2,1]*py[:,None] + bc[2]
        bounds = (px.min(), px.max(), py.min(), py.max())
        self.vector_plane = VectorPlane(vectors, self.box_vectors[2],
                                    bc, bounds, self.image,
                                    self.source._x_vec, self.source._y_vec)
        self.vp_bounds = bounds
        self.vectors = vectors

data_object_registry["volume_rendering"] = VolumeRendering

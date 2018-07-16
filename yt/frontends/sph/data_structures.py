"""
Data structures for SPH frontends.




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import os

from yt.data_objects.static_output import \
    ParticleDataset
from yt.funcs import mylog
from yt.geometry.particle_geometry_handler import \
    ParticleIndex

class SPHDataset(ParticleDataset):
    default_kernel_name = "cubic"
    _sph_smoothing_styles = ["scatter", "gather"]
    _sph_smoothing_style = "scatter"
    _num_neighbors = 32

    def __init__(self, filename, dataset_type=None, file_style=None,
                 units_override=None, unit_system="cgs",
                 index_order=None, index_filename=None,
                 kdtree_filename=None, kernel_name=None):
        if kernel_name is None:
            self.kernel_name = self.default_kernel_name
        else:
            self.kernel_name = kernel_name
        self.kdtree_filename = kdtree_filename
        super(SPHDataset, self).__init__(
            filename, dataset_type=dataset_type, file_style=file_style,
            units_override=units_override, unit_system=unit_system,
            index_order=index_order, index_filename=index_filename)

    @property
    def num_neighbors(self):
        return self._num_neighbors

    @num_neighbors.setter
    def num_neighbors(self, value):
        if value < 0:
            raise ValueError("Negative value not allowed: %s" % value)
        self._num_neighbors = value

    @property
    def sph_smoothing_style(self):
        return self._sph_smoothing_style

    @sph_smoothing_style.setter
    def sph_smoothing_style(self, value):
        if value not in self._sph_smoothing_styles:
            raise ValueError("Smoothing style not implemented: %s, please "
                             "select one of the following: " % value,
                             self._sph_smoothing_styles)

        self._sph_smoothing_style = value

    def add_smoothed_particle_field(self, smooth_field,
                                    method="volume_weighted", nneighbors=64,
                                    kernel_name=None):
        """Add a new smoothed particle field

        Creates a new smoothed field based on the particle *smooth_field*.

        Parameters
        ----------

        smooth_field : tuple
           The field name tuple of the particle field the smoothed field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        method : string, default 'volume_weighted'
           The particle smoothing method to use. Can only be 'volume_weighted'
           for now.
        nneighbors : int, default 64
            The number of neighbors to examine during the process.
        kernel_name : string or None, default None
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`. If left as None,
            :attr:`~yt.frontends.sph.data_structures.SPHDataset.kernel_name`
            will be used.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        if kernel_name is None:
            kernel_name = self.kernel_name
        return super(SPHDataset, self).add_smoothed_particle_field(
            smooth_field=smooth_field, method=method, nneighbors=nneighbors,
            kernel_name=kernel_name
        )

class SPHParticleIndex(ParticleIndex):
    def _initialize_index(self):
        ds = self.dataset

        ds._file_hash = self._generate_hash()

        if hasattr(self.io, '_generate_smoothing_length'):
            self.io._generate_smoothing_length(self.data_files, self.kdtree)

        super(SPHParticleIndex, self)._initialize_index()

    def _generate_kdtree(self, fname):
        from cykdtree import PyKDTree
        if os.path.exists(fname):
            mylog.info('Loading KDTree from %s' % os.path.basename(fname))
            kdtree = PyKDTree.from_file(fname)
            if kdtree.data_version != self.ds._file_hash:
                mylog.info('Detected hash mismatch, regenerating KDTree')
            else:
                self._kdtree = kdtree
                return
        positions = []
        for data_file in self.data_files:
            for _, ppos in self.io._yield_coordinates(
                    data_file, needed_ptype=self.ds._sph_ptype):
                positions.append(ppos)
        if positions == []:
            self._kdtree = None
            return
        positions = np.concatenate(positions)
        mylog.info('Allocating KDTree for %s particles' % positions.shape[0])
        self._kdtree = PyKDTree(
            positions.astype('float64'),
            left_edge=self.ds.domain_left_edge,
            right_edge=self.ds.domain_right_edge,
            periodic=np.array(self.ds.periodicity),
            leafsize=2*int(self.ds._num_neighbors),
            data_version=self.ds._file_hash
        )
        if fname is not None:
            self._kdtree.save(fname)

    @property
    def kdtree(self):
        if hasattr(self, '_kdtree'):
            return self._kdtree

        ds = self.ds

        if getattr(ds, 'kdtree_filename', None) is None:
            if os.path.exists(ds.parameter_filename):
                fname = ds.parameter_filename + ".kdtree"
            else:
                # we don't want to write to disk for in-memory data
                fname = None
        else:
            fname = ds.kdtree_filename

        self._generate_kdtree(fname)

        return self._kdtree

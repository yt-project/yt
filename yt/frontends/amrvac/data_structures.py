"""
AMRVAC data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import stat
import numpy as np
import weakref

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from .fields import AMRVACFieldInfo
from .misc import AMRVACDatReader


class AMRVACGrid(AMRGridPatch):
    """devnote : a patch represent part of a block. The hierarchy/index is a collection of patches"""
    _id_offset = 0

    def __init__(self, id, index, level):
        super(AMRVACGrid, self).__init__(
            id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "AMRVACGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class AMRVACHierarchy(GridIndex):
    grid = AMRVACGrid

    def __init__(self, ds, dataset_type="amrvac"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64

        # init everything to make it clear what's in there
        self.field_list = []
        self.num_grids = None

        super(AMRVACHierarchy, self).__init__(ds, dataset_type)


    def _detect_output_fields(self):
        # devnote: probably should distinguish gas and dust fields here
        # through the "field type" tag, which here is using self.dataset_type
        self.field_list = [(self.dataset_type, f) for f in self.dataset.parameters["w_names"]]

    def _count_grids(self):
        self.num_grids = self.dataset.parameters["nleafs"]

    def _parse_index(self):
        with open(self.dataset.parameter_filename, 'rb') as df:
            #devnote: here I'm loading everything in the RAM, defeating the purpose
            # this is a tmp workaround
            leaves_dat = AMRVACDatReader.get_block_data(df)
        header = self.dataset.parameters

        #all of these are (ndim) arrays
        domain_width = header['xmax'] - header['xmin']
        nblocks = header['domain_nx'] / header['block_nx']
        block_width = domain_width / nblocks

        self.grids = np.empty(len(leaves_dat), dtype='object')

        # those are YTarray instances, already initialized with proper shape
        self.grid_left_edge[:]  = 0.0
        self.grid_right_edge[:] = 1.0

        for ip, patch in enumerate(leaves_dat):
            patch_width = block_width / 2**(patch['lvl']-1)
            patch_left_edge  = header['xmin'] + (patch['ix']-1) / 2**(patch['lvl']) * domain_width
            patch_right_edge = header['xmin'] + (patch['ix'])   / 2**(patch['lvl']) * domain_width
            for idim, ledge in enumerate(patch_left_edge):
                # workaround the variable dimensionality of input data
                # missing values are left to init values (0 ?)
                self.grid_left_edge[ip,idim]  = patch_left_edge[idim]
                self.grid_right_edge[ip,idim] = patch_right_edge[idim]
                self.grid_dimensions[ip,idim] = patch['w'].shape[idim]
            self.grids[ip] = self.grid(id=ip, index=self, level=patch['lvl'])

        levels = np.array([leaf["lvl"] for leaf in leaves_dat])
        self.grid_levels = levels.reshape(self.num_grids, 1)
        self.max_level = self.dataset.parameters["levmax"]
        assert self.dataset.parameters["levmax"] == max(levels)

    def _populate_grid_objects(self):
        for g in self.grids:
            # set up Children and Parent...
            pass
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()



class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(self, filename, dataset_type='amrvac',
                 storage_filename=None,
                 units_override=None):
        self.fluid_types += ('amrvac',) #devnote: input 'gas', 'dust' here ?
        super(AMRVACDataset, self).__init__(filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename
        # refinement factor between a grid and its subgrid
        # self.refine_by = 2

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        #devnote: I'm using the default code because there need to be something but
        # this needs revising
        self.length_unit = self.quan(1.0, "cm")
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")

    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        #   self.parameters             <= full of code-specific items of use
        #   self.domain_left_edge       <= array of float64                         OK
        #   self.domain_right_edge      <= array of float64                         OK
        #   self.dimensionality         <= int                                      OK
        #   self.domain_dimensions      <= array of int64                           TODO
        #   self.periodicity            <= three-element tuple of booleans          TODO
        #   self.current_time           <= simulation time in code units            OK
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1                          OK
        #   self.current_redshift           <= float                                OK
        #   self.omega_lambda               <= float                                OK
        #   self.omega_matter               <= float                                OK
        #   self.hubble_constant            <= float                                OK
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        with open(self.parameter_filename, 'rb') as df:
            self.parameters = AMRVACDatReader.get_header(df)

        self.current_time   = self.parameters['time']
        self.dimensionality = self.parameters['ndim'] #devnote, warining : ndir != ndim
        #self.domain_dimensions = self.parameters[''].astype('int64') #devnote: number of cells, or grids ??

        dle = np.zeros(3)
        dre = np.ones(3)
        for idim in range(self.dimensionality):
            dle[idim] = self.parameters['xmin'][idim]
            dre[idim] = self.parameters['xmax'][idim]

        self.domain_left_edge = dle
        self.domain_right_edge = dre

        #devnote: these could be made optional if needed
        self.cosmological_simulation = 0
        self.current_redshift        = 0.0
        self.omega_matter            = 0.0
        self.omega_lambda            = 0.0
        self.hubble_constant         = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        validation = False
        try:
            with open(args[0], 'rb') as fi:
                assert 'rho' in fi.readline().decode('latin-1')
            validation = True
        finally:
            return validation

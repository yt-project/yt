"""
FLASH-specific data structures



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
from yt.data_objects.static_output import \
    Dataset, ParticleFile
from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.geometry.particle_geometry_handler import \
    ParticleIndex
from yt.utilities.file_handler import \
    HDF5FileHandler
from yt.utilities.physical_ratios import cm_per_mpc
from .fields import FLASHFieldInfo

class FLASHGrid(AMRGridPatch):
    _id_offset = 1
    #__slots__ = ["_level_id", "stop_index"]
    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename = index.index_filename,
                              index = index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "FLASHGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

class FLASHHierarchy(GridIndex):

    grid = FLASHGrid
    _preload_implemented = True
    
    def __init__(self,ds,dataset_type='flash_hdf5'):
        self.dataset_type = dataset_type
        self.field_indexes = {}
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self._handle = ds._handle
        self._particle_handle = ds._particle_handle
        self.float_type = np.float64
        GridIndex.__init__(self,ds,dataset_type)

    def _initialize_data_storage(self):
        pass

    def _detect_output_fields(self):
        self.field_list = [("flash", s.decode("ascii","ignore"))
                           for s in self._handle["/unknown names"][:].flat]
        if ("/particle names" in self._particle_handle):
            self.field_list += [("io", "particle_" +
                                    s[0].decode("ascii","ignore").strip())
                                for s in self._particle_handle["/particle names"][:]]
    
    def _count_grids(self):
        try:
            self.num_grids = self.dataset._find_parameter(
                "integer", "globalnumblocks", True)
        except KeyError:
            try:
                self.num_grids = \
                    self._handle['simulation parameters']['total blocks'][0]
            except KeyError:
                self.num_grids = self._handle["/simulation parameters"][0][0]
        
    def _parse_index(self):
        f = self._handle # shortcut
        ds = self.dataset # shortcut
        f_part = self._particle_handle # shortcut
        
        # Initialize to the domain left / domain right
        ND = self.dataset.dimensionality
        DLE = self.dataset.domain_left_edge
        DRE = self.dataset.domain_right_edge
        for i in range(3):
            self.grid_left_edge[:,i] = DLE[i]
            self.grid_right_edge[:,i] = DRE[i]
        # We only go up to ND for 2D datasets
        self.grid_left_edge[:,:ND] = f["/bounding box"][:,:ND,0]
        self.grid_right_edge[:,:ND] = f["/bounding box"][:,:ND,1]
        # Move this to the parameter file
        try:
            nxb = ds.parameters['nxb']
            nyb = ds.parameters['nyb']
            nzb = ds.parameters['nzb']
        except KeyError:
            nxb, nyb, nzb = [int(f["/simulation parameters"]['n%sb' % ax])
                              for ax in 'xyz']
        self.grid_dimensions[:] *= (nxb, nyb, nzb)
        try:
            self.grid_particle_count[:] = f_part["/localnp"][:][:,None]
        except KeyError:
            self.grid_particle_count[:] = 0.0
        self._particle_indices = np.zeros(self.num_grids + 1, dtype='int64')
        if self.num_grids > 1:
            np.add.accumulate(self.grid_particle_count.squeeze(),
                              out=self._particle_indices[1:])
        else:
            self._particle_indices[1] = self.grid_particle_count.squeeze()
        # This will become redundant, as _prepare_grid will reset it to its
        # current value.  Note that FLASH uses 1-based indexing for refinement
        # levels, but we do not, so we reduce the level by 1.
        self.grid_levels.flat[:] = f["/refine level"][:][:] - 1
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i+1, self, self.grid_levels[i,0])
        

        # This is a possibly slow and verbose fix, and should be re-examined!
        rdx = (self.dataset.domain_width /
                self.dataset.domain_dimensions)
        nlevels = self.grid_levels.max()
        dxs = np.ones((nlevels+1,3),dtype='float64')
        for i in range(nlevels+1):
            dxs[i,:ND] = rdx[:ND]/self.dataset.refine_by**i
       
        if ND < 3:
            dxs[:,ND:] = rdx[ND:]

        # Because we don't care about units, we're going to operate on views.
        gle = self.grid_left_edge.ndarray_view()
        gre = self.grid_right_edge.ndarray_view()
        geom = self.dataset.geometry
        if geom != 'cartesian' and ND < 3:
            if geom == 'spherical' and ND < 2:
                gle[:,1] = 0.0
                gre[:,1] = np.pi
            gle[:,2] = 0.0
            gre[:,2] = 2.0 * np.pi
            return

        # Now, for cartesian data.
        for i in range(self.num_grids):
            dx = dxs[self.grid_levels[i],:]
            gle[i][:ND] = np.rint(gle[i][:ND]/dx[0][:ND])*dx[0][:ND]
            gre[i][:ND] = np.rint(gre[i][:ND]/dx[0][:ND])*dx[0][:ND]

    def _populate_grid_objects(self):
        ii = np.argsort(self.grid_levels.flat)
        gid = self._handle["/gid"][:]
        first_ind = -(self.dataset.refine_by**self.dataset.dimensionality)
        for g in self.grids[ii].flat:
            gi = g.id - g._id_offset
            # FLASH uses 1-indexed group info
            g.Children = [self.grids[i - 1] for i in gid[gi,first_ind:] if i > -1]
            for g1 in g.Children:
                g1.Parent = g
            g._prepare_grid()
            g._setup_dx()
        if self.dataset.dimensionality < 3:
            DD = (self.dataset.domain_right_edge[2] -
                  self.dataset.domain_left_edge[2])
            for g in self.grids:
                g.dds[2] = DD
        if self.dataset.dimensionality < 2:
            DD = (self.dataset.domain_right_edge[1] -
                  self.dataset.domain_left_edge[1])
            for g in self.grids:
                g.dds[1] = DD
        self.max_level = self.grid_levels.max()

class FLASHDataset(Dataset):
    _index_class = FLASHHierarchy
    _field_info_class = FLASHFieldInfo
    _handle = None
    
    def __init__(self, filename, dataset_type='flash_hdf5',
                 storage_filename = None,
                 particle_filename = None, 
                 units_override = None,
                 unit_system = "cgs"):

        self.fluid_types += ("flash",)
        if self._handle is not None: return
        self._handle = HDF5FileHandler(filename)

        self.particle_filename = particle_filename

        if self.particle_filename is None:
            # try to guess the particle filename
            try:
                self._particle_handle = HDF5FileHandler(filename.replace('plt_cnt', 'part'))
                self.particle_filename = filename.replace('plt_cnt', 'part')
                mylog.info('Particle file found: %s' % self.particle_filename.split('/')[-1])
            except IOError:
                self._particle_handle = self._handle
        else:
            # particle_filename is specified by user
            try:
                self._particle_handle = HDF5FileHandler(self.particle_filename)
            except:
                raise IOError(self.particle_filename)
        # Check if the particle file has the same time
        if self._particle_handle != self._handle:
            part_time = self._particle_handle.handle.get('real scalars')[0][1]
            plot_time = self._handle.handle.get('real scalars')[0][1]
            if not np.isclose(part_time, plot_time):
                self._particle_handle = self._handle
                mylog.warning('%s and %s are not at the same time. ' % (self.particle_filename, filename) +
                              'This particle file will not be used.')

        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.refine_by = 2

        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename

        self.parameters["HydroMethod"] = 'flash' # always PPM DE
        self.parameters["Time"] = 1. # default unit is 1...
        
    def _set_code_unit_attributes(self):

        if 'unitsystem' in self.parameters:
            if self['unitsystem'].lower() == "cgs":
                b_factor = 1.0
            elif self['unitsystem'].lower() == "si":
                b_factor = np.sqrt(4*np.pi/1e7)
            elif self['unitsystem'].lower() == "none":
                b_factor = np.sqrt(4*np.pi)
            else:
                raise RuntimeError("Runtime parameter unitsystem with "
                                   "value %s is unrecognized" % self['unitsystem'])
        else:
            b_factor = 1.
        if self.cosmological_simulation == 1:
            length_factor = 1.0 / (1.0 + self.current_redshift)
            temperature_factor = 1.0 / (1.0 + self.current_redshift)**2
        else:
            length_factor = 1.0
            temperature_factor = 1.0

        setdefaultattr(self, 'magnetic_unit', self.quan(b_factor, "gauss"))
        setdefaultattr(self, 'length_unit', self.quan(length_factor, "cm"))
        setdefaultattr(self, 'mass_unit', self.quan(1.0, "g"))
        setdefaultattr(self, 'time_unit', self.quan(1.0, "s"))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, "cm/s"))
        setdefaultattr(
            self, 'temperature_unit', self.quan(temperature_factor, "K"))

    def set_code_units(self):
        super(FLASHDataset, self).set_code_units()

    def _find_parameter(self, ptype, pname, scalar = False):
        nn = "/%s %s" % (ptype,
                {False: "runtime parameters", True: "scalars"}[scalar])
        if nn not in self._handle: raise KeyError(nn)
        for tpname, pval in zip(self._handle[nn][:,'name'],
                                self._handle[nn][:,'value']):
            if tpname.decode("ascii","ignore").strip() == pname:
                if hasattr(pval, "decode"):
                    pval = pval.decode("ascii", "ignore")
                if ptype == "string":
                    return pval.strip()
                else:
                    return pval
        raise KeyError(pname)

    def _parse_parameter_file(self):
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        if "file format version" in self._handle:
            self._flash_version = int(
                self._handle["file format version"][:])
        elif "sim info" in self._handle:
            self._flash_version = int(
                self._handle["sim info"][:]["file format version"])
        else:
            raise RuntimeError("Can't figure out FLASH file version.")
        # First we load all of the parameters
        hns = ["simulation parameters"]
        # note the ordering here is important: runtime parameters should
        # ovewrite scalars with the same name.
        for ptype in ['scalars', 'runtime parameters']:
            for vtype in ['integer', 'real', 'logical', 'string']:
                hns.append("%s %s" % (vtype, ptype))
        if self._flash_version > 7:
            for hn in hns:
                if hn not in self._handle:
                    continue
                for varname, val in zip(self._handle[hn][:,'name'],
                                        self._handle[hn][:,'value']):
                    vn = varname.strip()
                    if hn.startswith("string") :
                        pval = val.strip()
                    else :
                        pval = val
                    if vn in self.parameters and self.parameters[vn] != pval:
                        mylog.info("{0} {1} overwrites a simulation "
                                   "scalar of the same name".format(hn[:-1],vn))
                    if hasattr(pval, 'decode'):
                        pval = pval.decode("ascii", "ignore")
                    self.parameters[vn.decode("ascii", "ignore")] = pval
        if self._flash_version == 7:
            for hn in hns:
                if hn not in self._handle:
                    continue
                if hn is 'simulation parameters':
                    zipover = ((name, self._handle[hn][name][0])
                               for name in self._handle[hn].dtype.names)
                else:
                    zipover = zip(self._handle[hn][:,'name'],self._handle[hn][:,'value'])
                for varname, val in zipover:
                    vn = varname.strip()
                    if hasattr(vn, 'decode'):
                        vn = vn.decode("ascii", "ignore")
                    if hn.startswith("string"):
                        pval = val.strip()
                    else:
                        pval = val
                    if vn in self.parameters and self.parameters[vn] != pval:
                        mylog.info("{0} {1} overwrites a simulation "
                                   "scalar of the same name".format(hn[:-1],vn))
                    if hasattr(pval, 'decode'):
                        pval = pval.decode("ascii", "ignore")
                    self.parameters[vn] = pval
        
        # Determine block size
        try:
            nxb = self.parameters["nxb"]
            nyb = self.parameters["nyb"]
            nzb = self.parameters["nzb"]
        except KeyError:
            nxb, nyb, nzb = [int(self._handle["/simulation parameters"]['n%sb' % ax])
                             for ax in 'xyz'] # FLASH2 only!
        
        # Determine dimensionality
        try:
            dimensionality = self.parameters["dimensionality"]
        except KeyError:
            dimensionality = 3
            if nzb == 1: dimensionality = 2
            if nyb == 1: dimensionality = 1
            if dimensionality < 3:
                mylog.warning("Guessing dimensionality as %s", dimensionality)
        
        self.dimensionality = dimensionality

        self.geometry = self.parameters["geometry"]
        # Determine base grid parameters
        if 'lrefine_min' in self.parameters.keys(): # PARAMESH
            nblockx = self.parameters["nblockx"]
            nblocky = self.parameters["nblocky"]
            nblockz = self.parameters["nblockz"]
        else: # Uniform Grid
            nblockx = self.parameters["iprocs"]
            nblocky = self.parameters["jprocs"]
            nblockz = self.parameters["kprocs"]

        # In case the user wasn't careful
        if self.dimensionality <= 2: nblockz = 1
        if self.dimensionality == 1: nblocky = 1

        # Determine domain boundaries
        self.domain_left_edge = np.array(
            [self.parameters["%smin" % ax] for ax in 'xyz']).astype("float64")
        self.domain_right_edge = np.array(
            [self.parameters["%smax" % ax] for ax in 'xyz']).astype("float64")
        if self.dimensionality < 3:
            for d in [dimensionality]+list(range(3-dimensionality)):
                if self.domain_left_edge[d] == self.domain_right_edge[d]:
                    mylog.warning('Identical domain left edge and right edges '
                                  'along dummy dimension (%i), attempting to read anyway' % d)
                    self.domain_right_edge[d] = self.domain_left_edge[d]+1.0
        if self.dimensionality < 3 and self.geometry == "cylindrical":
            mylog.warning("Extending theta dimension to 2PI + left edge.")
            self.domain_right_edge[2] = self.domain_left_edge[2] + 2*np.pi
        elif self.dimensionality < 3 and self.geometry == "polar":
            mylog.warning("Extending theta dimension to 2PI + left edge.")
            self.domain_right_edge[1] = self.domain_left_edge[1] + 2*np.pi
        elif self.dimensionality < 3 and self.geometry == "spherical":
            mylog.warning("Extending phi dimension to 2PI + left edge.")
            self.domain_right_edge[2] = self.domain_left_edge[2] + 2*np.pi
        if self.dimensionality == 1 and self.geometry == "spherical":
            mylog.warning("Extending theta dimension to PI + left edge.")
            self.domain_right_edge[1] = self.domain_left_edge[1] + np.pi
        self.domain_dimensions = \
            np.array([nblockx*nxb,nblocky*nyb,nblockz*nzb])

        # Try to determine Gamma
        try:
            self.gamma = self.parameters["gamma"]
        except:
            mylog.info("Cannot find Gamma")
            pass

        # Get the simulation time
        self.current_time = self.parameters["time"]

        # Determine if this is a periodic box
        p = [self.parameters.get("%sl_boundary_type" % ax, None) == "periodic" for ax in 'xyz']
        self.periodicity = tuple(p)

        # Determine cosmological parameters.
        try: 
            self.parameters["usecosmology"]
            self.cosmological_simulation = 1
            self.current_redshift = 1.0/self.parameters['scalefactor'] - 1.0
            self.omega_lambda = self.parameters['cosmologicalconstant']
            self.omega_matter = self.parameters['omegamatter']
            self.hubble_constant = self.parameters['hubbleconstant']
            self.hubble_constant *= cm_per_mpc * 1.0e-5 * 1.0e-2 # convert to 'h'
        except:
            self.current_redshift = self.omega_lambda = self.omega_matter = \
                self.hubble_constant = self.cosmological_simulation = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = HDF5FileHandler(args[0])
            if "bounding box" in fileh["/"].keys():
                return True
        except:
            pass
        return False

    @classmethod
    def _guess_candidates(cls, base, directories, files):
        candidates = [_ for _ in files if
                      ("_hdf5_plt_cnt_" in _) or
                      ("_hdf5_chk_" in _)]
        # Typically, Flash won't have nested outputs.
        return candidates, (len(candidates) == 0)

    def close(self):
        self._handle.close()

class FLASHParticleFile(ParticleFile):
    pass

class FLASHParticleDataset(FLASHDataset):
    _index_class = ParticleIndex
    over_refine_factor = 1
    filter_bbox = False
    _file_class = FLASHParticleFile

    def __init__(self, filename, dataset_type='flash_particle_hdf5',
                 storage_filename = None,
                 units_override = None,
                 n_ref = 64, unit_system = "cgs"):

        if self._handle is not None: return
        self._handle = HDF5FileHandler(filename)
        self.n_ref = n_ref
        self.refine_by = 2
        Dataset.__init__(self, filename, dataset_type, units_override=units_override,
                         unit_system=unit_system)
        self.storage_filename = storage_filename

    def _parse_parameter_file(self):
        # Let the superclass do all the work but then
        # fix the domain dimensions
        super(FLASHParticleDataset, self)._parse_parameter_file()
        nz = 1 << self.over_refine_factor
        domain_dimensions = np.zeros(3, "int32")
        domain_dimensions[:self.dimensionality] = nz
        self.domain_dimensions = domain_dimensions
        self.filename_template = self.parameter_filename
        self.file_count = 1

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = HDF5FileHandler(args[0])
            if "bounding box" not in fileh["/"].keys() \
                and "localnp" in fileh["/"].keys():
                return True
        except IOError:
            pass
        return False

    @classmethod
    def _guess_candidates(cls, base, directories, files):
        candidates = [_ for _ in files if "_hdf5_part_" in _]
        # Typically, Flash won't have nested outputs.
        return candidates, (len(candidates) == 0)

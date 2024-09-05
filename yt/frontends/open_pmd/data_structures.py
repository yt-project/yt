from os import listdir, path
from re import match
from typing import Optional

import numpy as np
from packaging.version import Version

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.data_objects.time_series import DatasetSeries
from yt.frontends.open_pmd.fields import OpenPMDFieldInfo
from yt.frontends.open_pmd.misc import get_component, \
    is_const_component, pad_to_threed
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.units.yt_array import YTArray
from yt.utilities.lib.misc_utilities import get_box_grids_level
from yt.utilities.file_handler import OpenPMDFileHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _openpmd_api as openpmd_api

ompd_known_versions = [Version(_) for _ in ("1.0.0", "1.0.1", "1.1.0")]
opmd_required_attributes = ["openPMD", "basePath"]


class OpenPMDGrid(AMRGridPatch):
    """Represents chunk of data on-disk.

    This defines the index and offset for every mesh and particle type.
    It also defines parents and children grids (still some problems here).
    """

    _id_offset = 0
    __slots__ = ["_level_id"]
    # Every particle species and mesh might have different indices and offsets

    ftypes: Optional[list[str]] = []
    ptypes: Optional[list[str]] = []
    findex = 0
    foffset = 0
    pindex = 0
    poffset = 0

    def __init__(self, gid, index, level = -1, fi = 0, fo = 0, pi = 0, po = None, ft = None, pt = None):
        AMRGridPatch.__init__(self, gid, filename=index.index_filename, index=index)
        if ft is None:
            ft = []
        if pt is None:
            pt = []
        #these indices and offset are now for on-disk purposes
        self.findex = fi
        self.foffset = fo
        self.pindex = pi
        self.poffset = po
        self.ftypes = ft
        self.ptypes = pt
        self._parent_id = []
        self._children_ids = []
        self.Level = level

    def __str__(self):
        return "OpenPMDGrid_%04i (%s)" % (self.id, self.ActiveDimensions)
    
    def _read_particles(self):
        """
        Return the number of particles for each species in this grid and the array
        of indices for particle sorting. This sorting is only used for datasets without particlePatches.
        
        """
        f = self.index.dataset._handle
        axes_labels = self.index.dataset._axes_labels
        #gle = self.LeftEdge.d[:len(axes_labels)][::-1]
        #gre = self.RightEdge.d[:len(axes_labels)][::-1]

        #quick check
        print(len(self.index.grids), len(f.particles[list(f.particles)[0]]['position'][axes_labels[0]].available_chunks()))

        assert len(self.index.grids) == len(f.particles[list(f.particles)[0]]['position'][axes_labels[0]].available_chunks())
        #for species in f.particles:
        #    
        #    part_dict[str(species)] = {'count': num_part, 'mask': np.where(index_array)[0]}

        """
        part_dict = {}
        for species in f.particles:
            if len(axes_labels) == 3:
                index_array = \
                        np.where(get_component(f.particles[species]['position'], axes_labels[0]) >= gle[0], True, False) *\
                        np.where(get_component(f.particles[species]['position'], axes_labels[0]) <  gre[0], True, False) *\
                        np.where(get_component(f.particles[species]['position'], axes_labels[1]) >= gle[1], True, False) *\
                        np.where(get_component(f.particles[species]['position'], axes_labels[1]) <  gre[1], True, False) *\
                        np.where(get_component(f.particles[species]['position'], axes_labels[2]) >= gle[2], True, False) *\
                        np.where(get_component(f.particles[species]['position'], axes_labels[2]) <  gre[2], True, False)
            elif len(axes_labels) == 2:
                index_array = \
                    np.where(get_component(f.particles[species]['position'], axes_labels[0]) >= gle[0], True, False) *\
                    np.where(get_component(f.particles[species]['position'], axes_labels[0]) <  gre[0], True, False) *\
                    np.where(get_component(f.particles[species]['position'], axes_labels[1]) >= gle[1], True, False) *\
                    np.where(get_component(f.particles[species]['position'], axes_labels[1]) <  gre[1], True, False)
            else:
                index_array = \
                    np.where(get_component(f.particles[species]['position'], axes_labels[0]) >= gle[0], True, False) *\
                    np.where(get_component(f.particles[species]['position'], axes_labels[0]) < gre[0], True, False)
            num_part = np.sum(index_array)
            part_dict[str(species)] = {'count': num_part, 'mask': np.where(index_array)[0]}
        """
        #self._pdata = part_dict

    def _setup_dx(self):
        # has already been read in and stored in index
        self.dds = self.index.ds.arr(self.index.level_dds[self.Level, :], "code_length")

    def _prepare_grid(self):
        super()._prepare_grid()
        h = self.index  # cache it
        my_ind = self.id - self._id_offset
        self.ActiveDimensions = h.grid_dimensions[my_ind]
        self.LeftEdge = h.grid_left_edge[my_ind]
        self.RightEdge = h.grid_right_edge[my_ind]

    @property
    def Parent(self):
        if len(self._parent_id) == 0:
            return None
        #self.Parent.dds
        return [self.index.grids[pid - self._id_offset] for pid in self._parent_id]
    
    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset] for cid in self._children_ids]
    
    @property
    def NumParticles(self):
        print(list(self.index.dataset._handle.particles['electrons']), 'PRINTED')
        if "particlesPatches" in self.index.dataset._handle.particles[\
                list(self.index.dataset._handle.particles)[0]]:
            for species in self.index.dataset._handle.particles:
                for patch in list(self.index.dataset._handle.particles[str(species)]
                                  ['particlePatches']):
                    mylog.warning("ParticlePatches is not yet supported as it has not yet been tested :(")
            print('PartPatches')
        else:
            self._read_particles()
        out_dict = {}
        for species_key in list(self._pdata.keys()):
            out_dict[species_key] = self._pdata[species_key]['count']
        return out_dict
    
class OpenPMDHierarchy(GridIndex):
    """Defines which fields and particles are created and read from disk.
    Furthermore it defines the characteristics of the grids.
    """

    grid = OpenPMDGrid

    def __init__(self, ds, dataset_type="openPMD"):
        self.dataset_type = dataset_type
        self.dataset = ds
        self.index_filename = ds.parameter_filename
        self.directory = path.dirname(self.index_filename)
        self.max_level = self.dataset._max_level
        GridIndex.__init__(self, ds, dataset_type)

    def _get_particle_type_counts(self):
        """Reads the active number of particles for every species.

        Returns
        -------
        dict
            keys are ptypes
            values are integer counts of the ptype
        """
        result = {}
        f = self.dataset._handle

        try:
            for ptype in self.ds.particle_types_raw:
                if str(ptype) == "io":
                    spec = list(f.particles)[0]
                else:
                    spec = ptype
                part = f.particles[spec]
                pos = part["position"][list(part["position"])[0]]
                if is_const_component(pos):
                    result[ptype] = pos.shape[0]  # relic, should there be a difference?
                else:
                    result[ptype] = pos.shape[0]
        except KeyError:
            result["io"] = 0

        return result

    def _detect_output_fields(self):
        """Populates ``self.field_list`` with native fields (mesh and particle) on disk.

        Each entry is a tuple of two strings. The first element is the on-disk fluid
        type or particle type. The second element is the name of the field in yt.
        This string is later used for accessing the data.
        Convention suggests that the on-disk fluid type should be "openPMD",
        the on-disk particle type (for a single species of particles) is "io"
        or (for multiple species of particles) the particle name on-disk.
        """
        f = self.dataset._handle
        mesh_fields = []
        try:
            # here we force higher level meshes to not appear on field list
            for mname in list(f.meshes)[:: self.max_level + 1]:
                try:
                    for axis in list(f.meshes[mname]):
                        mesh_fields.append(mname.replace("_", "-") + "_" + axis)
                except AttributeError:  # not sure if this would ever happen
                    # (i.e. scalar or constant component with no axes)
                    mesh_fields.append(mname.replace("_", "-"))
        except (KeyError, TypeError, AttributeError):
            pass
        self.field_list = [("openPMD", str(field)) for field in mesh_fields]

        particle_fields = []
        try:
            for pname in list(f.particles):
                species = list(f.particles[pname])
                for recname in species:
                    record = f.particles[pname][recname]
                    if is_const_component(record):
                        # Record itself (e.g. particle_mass) is constant
                        particle_fields.append(
                            pname.replace("_", "-") + "_" + recname.replace("_", "-")
                        )
                    elif "particlePatches" not in recname:
                        try:
                            # Create a field for every axis (x,y,z) of every
                            # property (position) of every species (electrons)
                            axes = list(record)
                            if str(recname) == "position":
                                recname = "positionCoarse"
                            if len(axes) > 1:
                                for axis in axes:
                                    particle_fields.append(
                                        pname.replace("_", "-")
                                        + "_"
                                        + recname.replace("_", "-")
                                        + "_"
                                        + axis
                                    )  # so we are doing all this, and then raising error which is why we get the noaxes field
                            else:
                                raise AttributeError  # in the case that is have no axes
                        except AttributeError:
                            # Record is a dataset, does not have axes (e.g. weighting) #electrons and ions momentum, positionCoarse, and positionOffset are here
                            particle_fields.append(
                                pname.replace("_", "-")
                                + "_"
                                + recname.replace("_", "-")
                            )
                            pass
                    else:
                        pass
            if len(list(f.particles)) > 1:
                # There is more than one particle species,
                # use the specific names as field types
                self.field_list.extend(
                    [
                        (
                            str(field).split("_")[0],
                            ("particle_" + "_".join(str(field).split("_")[1:])),
                        )
                        for field in particle_fields
                    ]
                )
            else:
                # Only one particle species, fall back to "io"
                self.field_list.extend(
                    [
                        ("io", ("particle_" + "_".join(str(field).split("_")[1:])))
                        for field in particle_fields
                    ]
                )
        except (KeyError, TypeError, AttributeError):
            pass

    def _count_grids(self):
        """Sets ``self.num_grids`` to be the total number of grids in the simulation.

        The number of grids is determined by their respective memory footprint.
        """
        f = self.dataset._handle
        self.meshshapes = {}
        self.numparts = {}
        self.num_grids = 0

        try:
            #here we assuming symmetry across meshes on a per_level basis
            for level, mname in enumerate(list(f.meshes)[: self.max_level + 1]):
                mesh = f.meshes[mname]
                if isinstance(mesh, openpmd_api.io.openpmd_api_cxx.Mesh):
                    if len(mesh[list(mesh)[0]].available_chunks()) > 0:
                        chunk_list = mesh[list(mesh)[0]].available_chunks()
                    else:
                        raise AttributeError
                    shape = tuple(mesh[list(mesh)[0]].shape) #per-level domain shape, not grid shape
                else:
                    raise AttributeError
                # currently we make one meshshape key/value dict pair per openpmdapi mesh
                spacing = tuple(mesh.grid_spacing)
                offset = tuple(mesh.grid_global_offset)
                unit_si = mesh.grid_unit_SI
                self.meshshapes[level] = (
                    chunk_list,
                    shape,
                    spacing,
                    offset,
                    unit_si,
                )
        except (TypeError, AttributeError):
            pass
        try:
            for pname in list(f.particles):
                species = f.particles[pname]
                if "particlePatches" in list(
                    species
                ):  # haven't been able to test this, likely doesn't work
                    for patch, size in enumerate(
                        species.particle_patches["numParticles"]
                    ):
                        self.numparts[f"{pname}#{patch}"] = size
                else:
                    axis = list(species["position"])[0]
                    if is_const_component(species["position"][axis]):
                        self.numparts[pname] = species["position"][axis].shape[0]
                    else:
                        self.numparts[pname] = species["position"][axis].shape[0]
                        # should these be the same? relic from old frontend
        except (KeyError, TypeError, AttributeError):
            pass

        # Limit values per grid by resulting memory footprint
        # TODO
        self.vpg = int(
            self.dataset.gridsize
        )  # 4Byte per value (f32) #havn't used this yet
        for chunk_ls, *_ in self.meshshapes.values():
            # ) #we are just using the amount of grids per level,
            # assuming equality between records of same level
            self.num_grids += len(chunk_ls)  # here we are not limiting by memory
        
        # Same goes for particle chunks if they are not inside particlePatches
        # what does this do?
        patches = {}
        no_patches = {}
        for k, v in self.numparts.items():
            if "#" in k:
                patches[k] = v
            else:
                no_patches[k] = v
        #Do we want empty grids for particles? I don't think so
        
        #for size in set(no_patches.values()):
        #    self.num_grids += int(np.ceil(size * self.vpg**-1))
        #for size in patches.values():
        #    self.num_grids += int(np.ceil(size * self.vpg**-1))

    def _parse_index(self):
        """Fills each grid with appropriate properties (extent, dimensions, ...)

        This calculates the properties of every OpenPMDGrid based on the total number of
        grids in the simulation. The domain is divided into ``self.num_grids`` based either on
        mesh_record_component.available_chunks() (only for ADIOS2 data) or memory contraints 
        (HDF5 data for now).

        Notes
        -----
        ``self.grid_dimensions`` is rounded to the nearest integer. Grid edges are
        calculated from this dimension. Grids with dimensions [0, 0, 0] are particle
        only. The others do not have any particles affiliated with them.
        """
        f = self.dataset._handle
        self.grid_levels.flat[:] = np.arange(self.max_level + 1)
        self.grids = np.empty(self.num_grids, dtype="object")

        grid_index_total = 0
        self.level_dds = np.empty((self.max_level + 1, 3)) #we always pad to three dimensions regardless of coordinate system
        # Mesh grids
        for level, mesh in self.meshshapes.items():
            (chunk_ls, shape, spacing, offset, unit_si) = mesh
            shape = np.asarray(shape)
            spacing = np.asarray(spacing)
            offset = np.asarray(offset)
            print('level', level, type(level))
            # Total dimension of this domain on a per-mesh-level basis!
            #domain_dimension = pad_to_threed(np.asarray(shape, dtype=np.int32), 1, self.dataset._geometry, self.dataset._axes_labels, self.dataset._data_order)

            for chunk in chunk_ls:  # convert chunks to grids!
                # dimension of individual chunks/grids
                #[::-1]?
                chunk_dim = pad_to_threed(np.array(chunk.extent, dtype = np.int32),
                                          1, self.dataset._geometry, self.dataset._axes_labels,  self.dataset._data_order)

                gle = (
                    np.array(chunk.offset) * unit_si * spacing + offset
                )  # to get to physical units, don't we need spacing?
                gre = (
                    np.array(chunk.offset) + np.array(chunk.extent)
                ) * unit_si * spacing + offset
                
                gle = pad_to_threed(gle, 0, self.dataset._geometry, self.dataset._axes_labels, self.dataset._data_order)
                gre = pad_to_threed(gre, 1, self.dataset._geometry, self.dataset._axes_labels,  self.dataset._data_order)
                # set things up

                self.grid_dimensions[grid_index_total] = chunk_dim
                self.grid_particle_count[grid_index_total] = 0 # does this change for particle patches?
                self.grid_left_edge[grid_index_total] = gle
                self.grid_right_edge[grid_index_total] = gre

                # we are hiding higher level access
                mesh_names = list(f.meshes)[0 :: self.max_level + 1]
                if len(list(f.particles))>1:
                    part_names = list(f.particles)
                else:
                    part_names = ["io"]

                chunk_index = np.array(chunk.offset, dtype=np.int32)
                chunk_extent = np.array(chunk.extent, dtype = np.int32)
                

                self.grids[grid_index_total] = self.grid(
                    grid_index_total,
                    self,
                    level,
                    fi = chunk_index,  # field index, openpmd-api chunk.offset
                    fo = chunk_extent,  # field offset, openpmd-api chunk.extent
                    ft = mesh_names,
                    pt = part_names,
                    #pi = 0, #what should these be, because particles might not be sorted?
                    #po = None,
                    #pd = num_parts #particle dictionary, mirroring _pdata dict from amrex frontend
                    
                )
                grid_index_total += 1
            print(unit_si)
            self.level_dds[level,:] = (gre-gle)/chunk_dim

        handled_ptypes = []

        # Particle grids
        for species, count in self.numparts.items():
            #HAS NOT BEEN TESTED
            if "#" in species:
                # This is a particlePatch, which should correspond with openPMD record component chunks
                """
                spec = species.split("#")
                patch = f.particles[spec[0]]["particlePatches"]  # don't know about this
                domain_dimension = np.ones(3, dtype=np.int32)
                for ind, axis in enumerate(list(patch["extent"])):
                    domain_dimension[ind] = patch["extent"][axis][int(spec[1])] #not sure about this
                num_grids = int(np.ceil(count * self.vpg**-1))
                gle = []
                for axis in patch["offset"]:
                    gle.append(
                        get_component(patch["offset"], axis, int(spec[1]), 1)[0]
                    )
                gle = np.asarray(gle)
                gle = np.append(gle, np.zeros(3 - len(gle)))
                gre = []
                for axis in patch["extent"]:
                    gre.append(
                        get_component(patch["extent"], axis, int(spec[1]), 1)[0]
                    )
                gre = np.asarray(gre)
                gre = np.append(gre, np.ones(3 - len(gre)))
                np.add(gle, gre, gre)
                npo = patch["numParticlesOffset"][()].item(int(spec[1]))
                particle_count = np.linspace(
                    npo, npo + count, num_grids + 1, dtype=np.int32
                )
                particle_names = [str(spec[0])]
                """

                """
                elif str(species) not in handled_ptypes:
                    domain_dimension = self.dataset.domain_dimensions
                    num_grids = int(np.ceil(count * self.vpg**-1))
                    gle = self.dataset.domain_left_edge
                    gre = self.dataset.domain_right_edge
                    particle_count = np.linspace(0, count, num_grids + 1, dtype=np.int32)
                    particle_names = []
                    for pname, size in self.numparts.items():
                        if size == count:
                            # Since this is not part of a particlePatch,
                            # we can include multiple same-sized ptypes
                            particle_names.append(str(pname))
                            handled_ptypes.append(str(pname))
                    """
            else:
                # A grid with this exact particle count has already been created
                continue
            for grid in np.arange(num_grids):
                """
                DONT MAKE NEW PARTICLE GRIDS, just add to existing field ones"""
                #self.grid_dimensions[grid_index_total] = domain_dimension
                #self.grid_left_edge[grid_index_total] = gle
                #self.grid_right_edge[grid_index_total] = gre

                self.grid_particle_count[grid_index_total] = (
                    particle_count[grid + 1] - particle_count[grid]
                ) * len(particle_names)
                self.grids[grid_index_total] = self.grid(
                    grid_index_total,
                    self,
                    0,
                    pi=particle_count[grid],
                    po=particle_count[grid + 1] - particle_count[grid],
                    pt=particle_names,
                )
                grid_index_total += 1
        print('ENDINDEX')
    
    def _populate_grid_objects(self):
        """This initializes all grids.
        Additionally, it should set up Children and Parent lists on each grid object.
        This almost works but calling 'self._reconstruct_parent_child()' before the loop 
        (with initializing the AMRGRidPatch baseclass) causes an attribute error in grid_patch.py:179
        AttributeError: 'OpenPMDGrid' object has no attribute 'LeftEdge'. 
        Calling 'self._reconstruct_parent_child()' after the loop causes an attribute error in grid_patch.py:68
        AttributeError: 'list' object has no attribute 'dds' because the Parent object is a list of grids and
        not a single grid as we allow overlapping edges for grids of different levels.

        """
        self.grids = np.array(self.grids, dtype="object")
        #self._reconstruct_parent_child()
        for i in np.arange(self.num_grids):
            self.grids[i]._prepare_grid()
            self.grids[i]._setup_dx()
        #self._reconstruct_parent_child()

    def _reconstruct_parent_child(self):
        if self.max_level == 0:
            return
        mask = np.empty(len(self.grids), dtype="int32")
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            print(i, 'index')
            get_box_grids_level(
                self.grid_left_edge[i, :],
                self.grid_right_edge[i, :],
                self.grid_levels[i].item() + 1,
                self.grid_left_edge,
                self.grid_right_edge,
                self.grid_levels,
                mask,
            )
            ids = np.where(mask.astype("bool"))  # where is a tuple
            grid._children_ids = ids[0] + grid._id_offset
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids):  # Second pass
            for child in grid.Children:
                child._parent_id.append(i + grid._id_offset)



class OpenPMDDataset(Dataset):
    """Contains all the required information of a single iteration of the simulation.

    Notes
    -----
    It is assumed that
    - all meshes cover the same region. Their resolution can be different.
    - all particles reside in this same region exclusively.
    - particle and mesh positions are *absolute* with respect to the simulation origin.
    """

    _load_requirements = ["openpmd_api"]
    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="openPMD",
        storage_filename=None,
        units_override=None,
        unit_system="mks",
        **kwargs,
    ):
        try:
            self._series_handle = OpenPMDFileHandler(filename)
            self._handle = self._series_handle.handle.iterations[
                list(self._series_handle.handle.iterations)[0]
            ]
        except TypeError:
            pass
        self.gridsize = kwargs.pop("open_pmd_virtual_gridsize", 10**9)
        self.standard_version = Version(self._series_handle.handle.openPMD)
        self.iteration = kwargs.pop("iteration", None)
        self._set_paths(self._series_handle, path.dirname(filename), self.iteration)
        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.storage_filename = storage_filename
        self.fluid_types += ("openPMD",)
        try:
            particles = tuple(str(c) for c in self._handle.particles)
            if len(particles) > 1:
                # Only use on-disk particle names if there is more than one species
                self.particle_types = particles
            mylog.debug("self.particle_types: %s", self.particle_types)
            self.particle_types_raw = self.particle_types
            self.particle_types = tuple(self.particle_types)
        except (KeyError, TypeError, AttributeError):
            pass

    def _set_paths(self, handle, path, iteration):
        """Parses relevant backend-paths out of ``handle``.
        Now this should be agnostic to backend format, just handling
        openpmd-api handle paths

        Parameters
        ----------
        handle : openpmd_api.openpmd_api_cxx.Series
        path : str
            (absolute) filepath for current hdf5, adios2, or (in the future) json container
        """
        iterations = []
        if iteration is None:
            iteration = list(handle.handle.iterations)[0]
            # moved this in
            encoding = str(handle.handle.iteration_encoding)
            if "Iteration_Encoding.group_based" in encoding:
                iterations = list(handle.handle.iterations)
                mylog.info(
                    "Single iteration found: %s", iteration
                )  # because of the groupbased dataset is_valid
            elif "Iteration_Encoding.file_based" in encoding:
                itformat = handle.handle.iteration_format
                regex = (
                    "^" + itformat.replace("%06T", "[0-9]+") + "$"
                )  # issues converting between yt and openpmd_api file patterns
                if path == "":
                    mylog.warning(
                        "For file based iterations, please use absolute file paths!"
                    )
                    pass
                for filename in listdir(path):
                    if match(regex, filename):
                        iterations.append(filename)
                mylog.info("Found %s iterations in directory", len(iterations))

            if len(iterations) == 0:
                mylog.warning("No iterations found!")
            if "Iteration_Encoding.group_based" in encoding and len(iterations) > 1:
                # this can happen only when a user doesn't use yt.load
                mylog.warning(
                    "More than one iteration found. Use OpenPMDDatasetSeries to view "
                    "other iterations. Loaded first iteration (%s)",
                    iteration,
                )
        # moved all this in
        self.base_path = f"/data/{iteration}/"
        try:  # we could find other ways to do this, but works to convert to iteration
            self.meshes_path = self._series_handle.handle.meshes_path
        except openpmd_api.io.ErrorNoSuchAttribute:
            if self.standard_version <= Version("1.1.0"):
                mylog.info(
                    "meshes_path not present in file. "
                    "Assuming file contains no meshes and has a domain extent of 1m^3!"
                )
                self.meshes_path = None
            else:
                raise
        try:
            self.particles_path = self._series_handle.handle.particles_path
        except openpmd_api.io.ErrorNoSuchAttribute:
            if self.standard_version <= Version("1.1.0"):
                mylog.info(
                    "particles_path not present in file."
                    " Assuming file contains no particles!"
                )
                self.particles_path = None
            else:
                raise

    def _set_code_unit_attributes(self):
        """Handle conversion between different physical units and the code units.

        Every dataset in openPMD can have different code <-> physical scaling.
        The individual factor is obtained by multiplying with "unitSI" reading getting
        data from disk.
        """
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "m/s"))
        setdefaultattr(self, "magnetic_unit", self.quan(1.0, "T"))

    def _parse_parameter_file(self):
        """Read in metadata describing the overall data on-disk."""
        f = self._handle  # this is an iteration
        self.parameters = 0
        self._periodicity = np.zeros(3, dtype="bool")
        self.refine_by = 1
        self.cosmological_simulation = 0
        max_lev = 0
        first_mesh = f.meshes[list(f.meshes)[0]]
        self._axes_labels = first_mesh.axis_labels
        self._data_order = first_mesh.data_order
        self._geometry = str(first_mesh.geometry)
        del(first_mesh)
        
        try:
            shapes = {}
            left_edges = {}
            right_edges = {}
            for mname in f.meshes:
                mesh = f.meshes[mname]
                #only access these three times, but they should all be the same?
                if isinstance(mesh, openpmd_api.io.openpmd_api_cxx.Mesh):
                    #all mesh axes should have same shape
                    shape = np.asarray(mesh[list(mesh)[0]].shape) 
                else:
                    shape = np.asarray(
                        mesh.shape
                    )  # if these aren't meshes, what are they?
                spacing = np.asarray(mesh.grid_spacing)
                offset = np.asarray(mesh.grid_global_offset)
                unit_si = np.asarray(mesh.grid_unit_SI)
                le = offset * unit_si
                if "lvl" in mname:
                    max_lev = max(max_lev, int(mname.split("lvl")[-1]))
                re = le + shape * unit_si * spacing
                shapes[mname] = shape #component_ordering(shape, self._geometry, self._data_order, self._axes_labels)
                left_edges[mname] = le #component_ordering(le, self._geometry, self._data_order, self._axes_labels)
                right_edges[mname] = re #component_ordering(re, self._geometry, self._data_order, self._axes_labels)
            # hidden as this is traditionally an index variable
            self._max_level = max_lev
            lowest_dim = np.min([len(i) for i in shapes.values()])
            #print([len(i) for i in shapes.values()])
            shapes = np.asarray([i[:lowest_dim] for i in shapes.values()])
            left_edges = np.asarray([i[:lowest_dim] for i in left_edges.values()])
            right_edges = np.asarray([i[:lowest_dim] for i in right_edges.values()])
            fs = []
            dle = []
            dre = []
            #Make these Row-ordered and then pad to three dimensions
            for i in np.arange(lowest_dim):
                fs.append(np.max(shapes.transpose()[i]).astype(int))
                dle.append(np.min(left_edges.transpose()[i]))
                dre.append(np.min(right_edges.transpose()[i]))
            self.dimensionality = len(fs)
            self.domain_dimensions = pad_to_threed(np.array(fs, dtype = np.int32), 1, self._geometry, self._axes_labels, self._data_order)
            self.domain_left_edge = pad_to_threed(np.array(dle, dtype = np.float64), 0, self._geometry, self._axes_labels, self._data_order)
            self.domain_right_edge = pad_to_threed(np.array(dre, dtype = np.float64), 1, self._geometry, self._axes_labels, self._data_order)
        except (KeyError, TypeError, AttributeError):
            #maybe we should throw an error for older versions of the standard.
            if self.standard_version <= Version("1.1.0"):
                self.dimensionality = 3
                self.domain_dimensions = np.ones(3, dtype=np.int32)
                self.domain_left_edge = np.zeros(3, dtype=np.float64)
                self.domain_right_edge = np.ones(3, dtype=np.float64)
            else:
                raise

        self.current_time = f.time * f.time_unit_SI

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        """Checks whether the supplied file can be read by this frontend."""
        if cls._missing_load_requirements():
            return False
        try:
            handle = openpmd_api.io.Series(
                filename, openpmd_api.io.Access_Type.read_only
            )
            for i in opmd_required_attributes:
                if i not in handle.attributes:
                    handle.close()
                    return False
            if Version(handle.openPMD) not in ompd_known_versions:
                handle.close()
                return False
            if "Iteration_Encoding.group_based" in str(handle.iteration_encoding):
                iteration = kwargs.pop("iteration", None)
                if len(list(handle.iterations)) > 1 and iteration is None:
                    handle.close()
                    return False
            handle.close()
            return True
        except (OSError, RuntimeError):
            return False


class OpenPMDDatasetSeries(DatasetSeries):
    _pre_outputs = ()
    _dataset_cls = OpenPMDDataset
    parallel = True
    setup_function = None
    mixed_dataset_types = False

    def __init__(self, filename):
        super().__init__([])
        # lets keep the yt patterns, then convert to %T and %06T for openpmd_api handling
        self.handle = openpmd_api.io.Series(
            filename, openpmd_api.io.Access_Type.read_only
        )
        if len(self.handle.iterations) < 2:
            mylog.info("Single iteration found, use OpenPMDDataset")
        self.filename = filename
        self._pre_outputs = sorted(
            np.asarray(list(self.handle.iterations), dtype="int64")
        )

    def __iter__(self):
        for it in self._pre_outputs:
            ds = self._load(it, **self.kwargs)
            self._setup_function(ds)
            yield ds

    def __getitem__(self, key):
        if isinstance(key, int):
            o = self._load(key)
            self._setup_function(o)
            return o
        else:
            raise KeyError(f"Unknown iteration {key}")

    def _load(self, it, **kwargs):
        return OpenPMDDataset(self.filename, iteration=it)


class OpenPMDGroupBasedDataset(Dataset):
    _load_requirements = ["openpmd_api"]
    _index_class = OpenPMDHierarchy
    _field_info_class = OpenPMDFieldInfo

    def __new__(cls, filename, *args, **kwargs):
        ret = object.__new__(OpenPMDDatasetSeries)
        ret.__init__(filename)
        mylog.info(
            "Group Based Dataset was cast to OpenPMDDatasetSeries.\n"
            "Inspect single interations with dataset_series[iteration]\n"
            "where iteration is in dataset_series._pre_outputs"
        )
        return ret

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if cls._missing_load_requirements():
            return False
        try:
            handle = openpmd_api.io.Series(
                filename, openpmd_api.io.Access_Type.read_only
            )
            for i in opmd_required_attributes:
                if i not in handle.attributes:
                    handle.close()
                    return False
            if Version(handle.openPMD) not in ompd_known_versions:
                handle.close()
                return False
            encoding = str(handle.iteration_encoding)
            iterations = list(handle.iterations)
            if "Iteration_Encoding.group_based" in encoding:
                iteration = kwargs.pop("iteration", None)
                # if iteration is given, not a dataset series
                if len(iterations) > 1 and iteration is None:
                    handle.close()
                    return True
            handle.close()
            return False
        except (OSError, RuntimeError):
            return False

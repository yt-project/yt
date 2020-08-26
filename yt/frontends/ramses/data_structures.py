import os
import weakref
from collections import defaultdict
from glob import glob

import numpy as np

from yt.arraytypes import blankRecordArray
from yt.data_objects.octree_subset import OctreeSubset
from yt.data_objects.particle_filters import add_particle_filter
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import RAMSESOctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.cython_fortran_utils import FortranFile as fpu
from yt.utilities.lib.cosmology_time import friedman
from yt.utilities.on_demand_imports import _f90nml as f90nml
from yt.utilities.physical_constants import kb, mp

from .definitions import field_aliases, particle_families, ramses_header
from .field_handlers import get_field_handlers
from .fields import _X, RAMSESFieldInfo
from .hilbert import get_cpu_list
from .io_utils import fill_hydro, read_amr
from .particle_handlers import get_particle_handlers


class RAMSESDomainFile:
    _last_mask = None
    _last_selector_id = None

    def __init__(self, ds, domain_id):
        self.ds = ds
        self.domain_id = domain_id

        num = os.path.basename(ds.parameter_filename).split(".")[0].split("_")[1]
        rootdir = ds.root_folder
        basedir = os.path.abspath(os.path.dirname(ds.parameter_filename))
        basename = "%s/%%s_%s.out%05i" % (basedir, num, domain_id)
        part_file_descriptor = f"{basedir}/part_file_descriptor.txt"
        if ds.num_groups > 0:
            igroup = ((domain_id - 1) // ds.group_size) + 1
            basename = "%s/group_%05i/%%s_%s.out%05i" % (
                rootdir,
                igroup,
                num,
                domain_id,
            )
        else:
            basename = "%s/%%s_%s.out%05i" % (basedir, num, domain_id)
        for t in ["grav", "amr"]:
            setattr(self, f"{t}_fn", basename % t)
        self._part_file_descriptor = part_file_descriptor
        self._read_amr_header()

        # Autodetect field files
        field_handlers = [FH(self) for FH in get_field_handlers() if FH.any_exist(ds)]
        self.field_handlers = field_handlers
        for fh in field_handlers:
            mylog.debug("Detected fluid type %s in domain_id=%s", fh.ftype, domain_id)
            fh.detect_fields(ds)
            # self._add_ftype(fh.ftype)

        # Autodetect particle files
        particle_handlers = [
            PH(ds, self) for PH in get_particle_handlers() if PH.any_exist(ds)
        ]
        self.particle_handlers = particle_handlers
        for ph in particle_handlers:
            mylog.debug(
                "Detected particle type %s in domain_id=%s", ph.ptype, domain_id
            )
            ph.read_header()
            # self._add_ptype(ph.ptype)

        # Load the AMR structure
        self._read_amr()

    _hydro_offset = None
    _level_count = None

    def __repr__(self):
        return "RAMSESDomainFile: %i" % self.domain_id

    @property
    def level_count(self):
        lvl_count = None
        for fh in self.field_handlers:
            fh.offset
            if lvl_count is None:
                lvl_count = fh.level_count.copy()
            else:
                lvl_count += fh._level_count
        return lvl_count

    @property
    def amr_file(self):
        if hasattr(self, "_amr_file"):
            self._amr_file.seek(0)
            return self._amr_file

        f = fpu(self.amr_fn)
        self._amr_file = f
        f.seek(0)
        return f

    def _read_amr_header(self):
        hvals = {}
        f = self.amr_file
        f.seek(0)

        for header in ramses_header(hvals):
            hvals.update(f.read_attrs(header))
        # For speedup, skip reading of 'headl' and 'taill'
        f.skip(2)
        hvals["numbl"] = f.read_vector("i")

        # That's the header, now we skip a few.
        hvals["numbl"] = np.array(hvals["numbl"]).reshape(
            (hvals["nlevelmax"], hvals["ncpu"])
        )
        f.skip()
        if hvals["nboundary"] > 0:
            f.skip(2)
            self.ngridbound = f.read_vector("i").astype("int64")
        else:
            self.ngridbound = np.zeros(hvals["nlevelmax"], dtype="int64")
        free_mem = f.read_attrs((("free_mem", 5, "i"),))  # NOQA
        ordering = f.read_vector("c")  # NOQA
        f.skip(4)
        # Now we're at the tree itself
        # Now we iterate over each level and each CPU.
        self.amr_header = hvals
        # update levelmax
        force_max_level, convention = self.ds._force_max_level
        if convention == "yt":
            force_max_level += self.ds.min_level + 1
        self.amr_header["nlevelmax"] = min(
            force_max_level, self.amr_header["nlevelmax"]
        )
        self.amr_offset = f.tell()
        self.local_oct_count = hvals["numbl"][
            self.ds.min_level :, self.domain_id - 1
        ].sum()
        self.total_oct_count = hvals["numbl"][self.ds.min_level :, :].sum(axis=0)

    def _read_amr(self):
        """Open the oct file, read in octs level-by-level.
           For each oct, only the position, index, level and domain
           are needed - its position in the octree is found automatically.
           The most important is finding all the information to feed
           oct_handler.add
        """
        self.oct_handler = RAMSESOctreeContainer(
            self.ds.domain_dimensions / 2,
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
        )
        root_nodes = self.amr_header["numbl"][self.ds.min_level, :].sum()
        self.oct_handler.allocate_domains(self.total_oct_count, root_nodes)
        mylog.debug(
            "Reading domain AMR % 4i (%0.3e, %0.3e)",
            self.domain_id,
            self.total_oct_count.sum(),
            self.ngridbound.sum(),
        )

        f = self.amr_file
        f.seek(self.amr_offset)

        min_level = self.ds.min_level
        max_level = read_amr(
            f, self.amr_header, self.ngridbound, min_level, self.oct_handler
        )

        self.max_level = max_level
        self.oct_handler.finalize()

        # Close AMR file
        f.close()

    def included(self, selector):
        if getattr(selector, "domain_id", None) is not None:
            return selector.domain_id == self.domain_id
        domain_ids = self.oct_handler.domain_identify(selector)
        return self.domain_id in domain_ids


class RAMSESDomainSubset(OctreeSubset):

    _domain_offset = 1
    _block_order = "F"

    _base_domain = None

    def __init__(
        self,
        base_region,
        domain,
        ds,
        over_refine_factor=1,
        num_ghost_zones=0,
        base_grid=None,
    ):
        super(RAMSESDomainSubset, self).__init__(
            base_region, domain, ds, over_refine_factor, num_ghost_zones
        )

        self._base_grid = base_grid

        if num_ghost_zones > 0:
            if not all(ds.periodicity):
                mylog.warning(
                    "Ghost zones will wrongly assume the domain to be periodic."
                )
            # Create a base domain *with no self._base_domain.fwidth
            base_domain = RAMSESDomainSubset(
                ds.all_data(), domain, ds, over_refine_factor
            )
            self._base_domain = base_domain
        elif num_ghost_zones < 0:
            raise RuntimeError(
                "Cannot initialize a domain subset with a negative number "
                "of ghost zones, was called with num_ghost_zones=%s" % num_ghost_zones
            )

    def _fill_no_ghostzones(self, fd, fields, selector, file_handler):
        ndim = self.ds.dimensionality
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        oct_handler = self.oct_handler
        all_fields = [f for ft, f in file_handler.field_list]
        fields = [f for ft, f in fields]
        data = {}
        cell_count = selector.count_oct_cells(self.oct_handler, self.domain_id)

        levels, cell_inds, file_inds = self.oct_handler.file_index_octs(
            selector, self.domain_id, cell_count
        )

        # Initializing data container
        for field in fields:
            data[field] = np.zeros(cell_count, "float64")

        cpu_list = [self.domain_id - 1]
        fill_hydro(
            fd,
            file_handler.offset,
            file_handler.level_count,
            cpu_list,
            levels,
            cell_inds,
            file_inds,
            ndim,
            all_fields,
            fields,
            data,
            oct_handler,
        )
        return data

    def _fill_with_ghostzones(
        self, fd, fields, selector, file_handler, num_ghost_zones
    ):
        ndim = self.ds.dimensionality
        ncpu = self.ds.parameters["ncpu"]
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        oct_handler = self.oct_handler
        all_fields = [f for ft, f in file_handler.field_list]
        fields = [f for ft, f in fields]
        tr = {}

        cell_count = (
            selector.count_octs(self.oct_handler, self.domain_id) * self.nz ** ndim
        )

        (
            levels,
            cell_inds,
            file_inds,
            domains,
        ) = self.oct_handler.file_index_octs_with_ghost_zones(
            selector, self.domain_id, cell_count
        )

        # Initializing data container
        for field in fields:
            tr[field] = np.zeros(cell_count, "float64")
        cpu_list = list(range(ncpu))
        fill_hydro(
            fd,
            file_handler.offset,
            file_handler.level_count,
            cpu_list,
            levels,
            cell_inds,
            file_inds,
            ndim,
            all_fields,
            fields,
            tr,
            oct_handler,
            domains=domains,
        )
        return tr

    @property
    def fwidth(self):
        fwidth = super(RAMSESDomainSubset, self).fwidth
        if self._num_ghost_zones > 0:
            fwidth = fwidth.reshape(-1, 8, 3)
            n_oct = fwidth.shape[0]
            # new_fwidth contains the fwidth of the oct+ghost zones
            # this is a constant array in each oct, so we simply copy
            # the oct value using numpy fancy-indexing
            new_fwidth = np.zeros((n_oct, self.nz ** 3, 3), dtype=fwidth.dtype)
            new_fwidth[:, :, :] = fwidth[:, 0:1, :]
            fwidth = new_fwidth.reshape(-1, 3)
        return fwidth

    @property
    def fcoords(self):
        num_ghost_zones = self._num_ghost_zones
        if num_ghost_zones == 0:
            return super(RAMSESDomainSubset, self).fcoords

        oh = self.oct_handler

        indices = oh.fill_index(self.selector).reshape(-1, 8)
        oct_inds, cell_inds = oh.fill_octcellindex_neighbours(self.selector)

        oct_inds = oct_inds.reshape(-1, 64)
        cell_inds = cell_inds.reshape(-1, 64)

        inds = indices[oct_inds, cell_inds]

        fcoords = self.ds.arr(oh.fcoords(self.selector)[inds].reshape(-1, 3), "unitary")

        return fcoords

    def fill(self, fd, fields, selector, file_handler):
        if self._num_ghost_zones == 0:
            return self._fill_no_ghostzones(fd, fields, selector, file_handler)
        else:
            return self._fill_with_ghostzones(
                fd, fields, selector, file_handler, self._num_ghost_zones
            )

    def retrieve_ghost_zones(self, ngz, fields, smoothed=False):
        if smoothed:
            mylog.warning(
                f"{self}.retrieve_ghost_zones was called with the "
                f"`smoothed` argument set to True. This is not supported, "
                "ignoring it."
            )
            smoothed = False

        try:
            new_subset = self._subset_with_gz
            mylog.debug(
                "Reusing previous subset with ghost zone for domain %s", self.domain_id
            )
        except AttributeError:
            new_subset = RAMSESDomainSubset(
                self.base_region,
                self.domain,
                self.ds,
                num_ghost_zones=ngz,
                base_grid=self,
            )
            self._subset_with_gz = new_subset

            # Cache the fields
            new_subset.get_data(fields)

        return new_subset


class RAMSESIndex(OctreeIndex):
    def __init__(self, ds, dataset_type="ramses"):
        self.fluid_field_list = ds._fields_in_file
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.max_level = None

        self.float_type = np.float64
        super(RAMSESIndex, self).__init__(ds, dataset_type)

    def _initialize_oct_handler(self):
        if self.ds._bbox is not None:
            cpu_list = get_cpu_list(self.dataset, self.dataset._bbox)
        else:
            cpu_list = range(self.dataset["ncpu"])

        self.domains = [RAMSESDomainFile(self.dataset, i + 1) for i in cpu_list]
        total_octs = sum(
            dom.local_oct_count for dom in self.domains  # + dom.ngridbound.sum()
        )
        force_max_level, convention = self.ds._force_max_level
        if convention == "yt":
            force_max_level += self.ds.min_level + 1
        self.max_level = min(
            force_max_level, max(dom.max_level for dom in self.domains)
        )
        self.num_grids = total_octs

    def _detect_output_fields(self):
        dsl = set([])

        # Get the detected particle fields
        for domain in self.domains:
            for ph in domain.particle_handlers:
                dsl.update(set(ph.field_offsets.keys()))

        self.particle_field_list = list(dsl)
        cosmo = self.ds.cosmological_simulation
        for f in self.particle_field_list[:]:
            if f[1] == "particle_birth_time" and cosmo:
                self.particle_field_list.append((f[0], "conformal_birth_time"))

        # Get the detected fields
        dsl = set([])
        for fh in self.domains[0].field_handlers:
            dsl.update(set(fh.field_list))
        self.fluid_field_list = list(dsl)

        self.field_list = self.particle_field_list + self.fluid_field_list

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            domains = [dom for dom in self.domains if dom.included(dobj.selector)]
            base_region = getattr(dobj, "base_region", dobj)
            if len(domains) > 1:
                mylog.debug("Identified %s intersecting domains", len(domains))
            subsets = [
                RAMSESDomainSubset(
                    base_region,
                    domain,
                    self.dataset,
                    num_ghost_zones=dobj._num_ghost_zones,
                )
                for domain in domains
            ]
            dobj._chunk_info = subsets
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g], None)

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None, cache=cache)

    def _initialize_level_stats(self):
        levels = sum([dom.level_count for dom in self.domains])
        desc = {"names": ["numcells", "level"], "formats": ["Int64"] * 2}
        max_level = self.dataset.min_level + self.dataset.max_level + 2
        self.level_stats = blankRecordArray(desc, max_level)
        self.level_stats["level"] = [i for i in range(max_level)]
        self.level_stats["numcells"] = [0 for i in range(max_level)]
        for level in range(self.dataset.min_level + 1):
            self.level_stats[level + 1]["numcells"] = 2 ** (
                level * self.dataset.dimensionality
            )
        for level in range(self.max_level + 1):
            self.level_stats[level + self.dataset.min_level + 1]["numcells"] = levels[
                level
            ]

    def _get_particle_type_counts(self):
        npart = 0
        npart = {k: 0 for k in self.ds.particle_types if k != "all"}
        for dom in self.domains:
            for fh in dom.particle_handlers:
                count = fh.local_particle_count
                npart[fh.ptype] += count

        return npart

    def print_stats(self):
        """
        Prints out (stdout) relevant information about the simulation

        This function prints information based on the fluid on the grids,
        and therefore does not work for DM only runs.
        """
        if not self.fluid_field_list:
            print("This function is not implemented for DM only runs")
            return

        self._initialize_level_stats()

        header = "%3s\t%14s\t%14s" % ("level", "# cells", "# cells^3")
        print(header)
        print(f"{len(header.expandtabs()) * '-'}")
        for level in range(self.dataset.min_level + self.dataset.max_level + 2):
            print(
                "% 3i\t% 14i\t% 14i"
                % (
                    level,
                    self.level_stats["numcells"][level],
                    np.ceil(self.level_stats["numcells"][level] ** (1.0 / 3)),
                )
            )
        print("-" * 46)
        print("   \t% 14i" % (self.level_stats["numcells"].sum()))
        print("\n")

        dx = self.get_smallest_dx()
        try:
            print("z = %0.8f" % (self.dataset.current_redshift))
        except Exception:
            pass
        print(
            "t = %0.8e = %0.8e s = %0.8e years"
            % (
                self.ds.current_time.in_units("code_time"),
                self.ds.current_time.in_units("s"),
                self.ds.current_time.in_units("yr"),
            )
        )
        print("\nSmallest Cell:")
        for item in ("Mpc", "pc", "AU", "cm"):
            print("\tWidth: %0.3e %s" % (dx.in_units(item), item))


class RAMSESDataset(Dataset):
    _index_class = RAMSESIndex
    _field_info_class = RAMSESFieldInfo
    gamma = 1.4  # This will get replaced on hydro_fn open

    def __init__(
        self,
        filename,
        dataset_type="ramses",
        fields=None,
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        extra_particle_fields=None,
        cosmological=None,
        bbox=None,
        max_level=None,
        max_level_convention=None,
    ):
        # Here we want to initiate a traceback, if the reader is not built.
        if isinstance(fields, str):
            fields = field_aliases[fields]
        """
        fields:
        An array of hydro variable fields in order of position in the
        hydro_XXXXX.outYYYYY file. If set to None, will try a default set of fields.

        extra_particle_fields:
        An array of extra particle variables in order of position in the
        particle_XXXXX.outYYYYY file.

        cosmological:
        If set to None, automatically detect cosmological simulation.
        If a boolean, force its value.
        """

        self._fields_in_file = fields
        # By default, extra fields have not triggered a warning
        self._warned_extra_fields = defaultdict(lambda: False)
        self._extra_particle_fields = extra_particle_fields
        self.force_cosmological = cosmological
        self._bbox = bbox

        self._force_max_level = self._sanitize_max_level(
            max_level, max_level_convention
        )

        # Infer if the output is organized in groups
        root_folder, group_folder = os.path.split(os.path.split(filename)[0])

        if group_folder == "group_00001":
            # Count the number of groups
            # note: we exclude the unlikely event that one of the group is actually a
            # file instad of a folder
            self.num_groups = len(
                [
                    _
                    for _ in glob(os.path.join(root_folder, "group_?????"))
                    if os.path.isdir(_)
                ]
            )
            self.root_folder = root_folder
        else:
            self.root_folder = os.path.split(filename)[0]
            self.num_groups = 0

        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

        # Add the particle types
        ptypes = []
        for PH in get_particle_handlers():
            if PH.any_exist(self):
                ptypes.append(PH.ptype)

        ptypes = tuple(ptypes)
        self.particle_types = self.particle_types_raw = ptypes

        # Add the fluid types
        for FH in get_field_handlers():
            FH.purge_detected_fields(self)
            if FH.any_exist(self):
                self.fluid_types += (FH.ftype,)

        self.storage_filename = storage_filename

    @staticmethod
    def _sanitize_max_level(max_level, max_level_convention):
        # NOTE: the initialisation of the dataset class sets
        #       self.min_level _and_ requires force_max_level
        #       to be set, so we cannot convert from to yt/ramses
        #       conventions
        if max_level is None and max_level_convention is None:
            return (2 ** 999, "yt")

        # Check max_level is a valid, positive integer
        if not isinstance(max_level, (int, np.integer)):
            raise TypeError(
                f"Expected `max_level` to be a positive integer, got {max_level} "
                f"with type {type(max_level)} instead."
            )
        if max_level < 0:
            raise ValueError(
                f"Expected `max_level` to be a positive integer, got {max_level} "
                "instead."
            )

        # Check max_level_convention is set and acceptable
        if max_level_convention is None:
            raise ValueError(
                f"Received `max_level`={max_level}, but no `max_level_convention`. "
                "Valid conventions are 'yt' and 'ramses'."
            )
        if max_level_convention not in ("ramses", "yt"):
            raise ValueError(
                f"Invalid convention {max_level_convention}. "
                "Valid choices are 'yt' and 'ramses'."
            )
        return (max_level, max_level_convention)

    def create_field_info(self, *args, **kwa):
        """Extend create_field_info to add the particles types."""
        super(RAMSESDataset, self).create_field_info(*args, **kwa)
        # Register particle filters
        if ("io", "particle_family") in self.field_list:
            for fname, value in particle_families.items():

                def loc(val):
                    def closure(pfilter, data):
                        filter = data[(pfilter.filtered_type, "particle_family")] == val
                        return filter

                    return closure

                add_particle_filter(
                    fname, loc(value), filtered_type="io", requires=["particle_family"]
                )

            for k in particle_families.keys():
                mylog.info("Adding particle_type: %s", k)
                self.add_particle_filter(f"{k}")

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]

    def _set_code_unit_attributes(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        # loading the units from the info file
        boxlen = self.parameters["boxlen"]
        length_unit = self.parameters["unit_l"]
        density_unit = self.parameters["unit_d"]
        time_unit = self.parameters["unit_t"]

        # calculating derived units (except velocity and temperature, done below)
        mass_unit = density_unit * length_unit ** 3
        magnetic_unit = np.sqrt(4 * np.pi * mass_unit / (time_unit ** 2 * length_unit))
        pressure_unit = density_unit * (length_unit / time_unit) ** 2

        # TODO:
        # Generalize the temperature field to account for ionization
        # For now assume an atomic ideal gas with cosmic abundances (x_H = 0.76)
        mean_molecular_weight_factor = _X ** -1

        setdefaultattr(self, "density_unit", self.quan(density_unit, "g/cm**3"))
        setdefaultattr(self, "magnetic_unit", self.quan(magnetic_unit, "gauss"))
        setdefaultattr(self, "pressure_unit", self.quan(pressure_unit, "dyne/cm**2"))
        setdefaultattr(self, "time_unit", self.quan(time_unit, "s"))
        setdefaultattr(self, "mass_unit", self.quan(mass_unit, "g"))
        setdefaultattr(
            self, "velocity_unit", self.quan(length_unit, "cm") / self.time_unit
        )
        temperature_unit = (
            self.velocity_unit ** 2 * mp * mean_molecular_weight_factor / kb
        )
        setdefaultattr(self, "temperature_unit", temperature_unit.in_units("K"))

        # Only the length unit get scales by a factor of boxlen
        setdefaultattr(self, "length_unit", self.quan(length_unit * boxlen, "cm"))

    def _parse_parameter_file(self):
        # hardcoded for now
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "ramses"
        self.parameters["Time"] = 1.0  # default unit is 1...

        # We now execute the same logic Oliver's code does
        rheader = {}

        def read_rhs(f, cast):
            line = f.readline().replace("\n", "")
            p, v = line.split("=")
            rheader[p.strip()] = cast(v.strip())

        with open(self.parameter_filename) as f:
            for _ in range(6):
                read_rhs(f, int)
            f.readline()
            for _ in range(11):
                read_rhs(f, float)
            f.readline()
            read_rhs(f, str)
            # This next line deserves some comment.  We specify a min_level that
            # corresponds to the minimum level in the RAMSES simulation.  RAMSES is
            # one-indexed, but it also does refer to the *oct* dimensions -- so
            # this means that a levelmin of 1 would have *1* oct in it.  So a
            # levelmin of 2 would have 8 octs at the root mesh level.
            self.min_level = rheader["levelmin"] - 1
            # Now we read the hilbert indices
            self.hilbert_indices = {}
            if rheader["ordering type"] == "hilbert":
                f.readline()  # header
                for _ in range(rheader["ncpu"]):
                    dom, mi, ma = f.readline().split()
                    self.hilbert_indices[int(dom)] = (float(mi), float(ma))

        if rheader["ordering type"] != "hilbert" and self._bbox is not None:
            raise NotImplementedError(
                "The ordering %s is not compatible with the `bbox` argument."
                % rheader["ordering type"]
            )
        self.parameters.update(rheader)
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_dimensions = np.ones(3, dtype="int32") * 2 ** (self.min_level + 1)
        self.domain_right_edge = np.ones(3, dtype="float64")
        # This is likely not true, but it's not clear
        # how to determine the boundary conditions
        self.periodicity = (True, True, True)

        if self.force_cosmological is not None:
            is_cosmological = self.force_cosmological
        else:
            # These conditions seem to always be true for non-cosmological datasets
            is_cosmological = not (
                rheader["time"] >= 0 and rheader["H0"] == 1 and rheader["aexp"] == 1
            )

        if not is_cosmological:
            self.cosmological_simulation = 0
            self.current_redshift = 0
            self.hubble_constant = 0
            self.omega_matter = 0
            self.omega_lambda = 0
        else:
            self.cosmological_simulation = 1
            self.current_redshift = (1.0 / rheader["aexp"]) - 1.0
            self.omega_lambda = rheader["omega_l"]
            self.omega_matter = rheader["omega_m"]
            self.hubble_constant = rheader["H0"] / 100.0  # This is H100

        force_max_level, convention = self._force_max_level
        if convention == "yt":
            force_max_level += self.min_level + 1
        self.max_level = min(force_max_level, rheader["levelmax"]) - self.min_level - 1

        if self.cosmological_simulation == 0:
            self.current_time = self.parameters["time"]
        else:
            self.tau_frw, self.t_frw, self.dtau, self.n_frw, self.time_tot = friedman(
                self.omega_matter,
                self.omega_lambda,
                1.0 - self.omega_matter - self.omega_lambda,
            )

            age = self.parameters["time"]
            iage = 1 + int(10.0 * age / self.dtau)
            iage = np.min([iage, self.n_frw // 2 + (iage - self.n_frw // 2) // 10])

            self.time_simu = self.t_frw[iage] * (age - self.tau_frw[iage - 1]) / (
                self.tau_frw[iage] - self.tau_frw[iage - 1]
            ) + self.t_frw[iage - 1] * (age - self.tau_frw[iage]) / (
                self.tau_frw[iage - 1] - self.tau_frw[iage]
            )

            self.current_time = (
                (self.time_tot + self.time_simu)
                / (self.hubble_constant * 1e7 / 3.08e24)
                / self.parameters["unit_t"]
            )

        if self.num_groups > 0:
            self.group_size = rheader["ncpu"] // self.num_groups

        # Read namelist.txt file (if any)
        self.read_namelist()

    def read_namelist(self):
        """Read the namelist.txt file in the output folder, if present"""
        namelist_file = os.path.join(self.root_folder, "namelist.txt")
        if os.path.exists(namelist_file):
            try:
                with open(namelist_file, "r") as f:
                    nml = f90nml.read(f)
            except ImportError as e:
                nml = f"An error occurred when reading the namelist: {str(e)}"
            except (ValueError, StopIteration) as e:
                mylog.warning(
                    "Could not parse `namelist.txt` file as it was malformed: %s", e
                )
                return

            self.parameters["namelist"] = nml

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not os.path.basename(args[0]).startswith("info_"):
            return False
        fn = args[0].replace("info_", "amr_").replace(".txt", ".out00001")
        return os.path.exists(fn)

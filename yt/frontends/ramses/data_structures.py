import os
import weakref
from collections import defaultdict
from functools import cached_property
from itertools import product
from pathlib import Path

import numpy as np

from yt.arraytypes import blankRecordArray
from yt.data_objects.index_subobjects.octree_subset import OctreeSubset
from yt.data_objects.particle_filters import add_particle_filter
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.oct_container import RAMSESOctreeContainer
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.cython_fortran_utils import FortranFile as fpu
from yt.utilities.lib.cosmology_time import t_frw, tau_frw
from yt.utilities.on_demand_imports import _f90nml as f90nml
from yt.utilities.physical_constants import kb, mp

from .definitions import (
    OUTPUT_DIR_EXP,
    OUTPUT_DIR_RE,
    STANDARD_FILE_RE,
    field_aliases,
    particle_families,
    ramses_header,
)
from .field_handlers import get_field_handlers
from .fields import RAMSESFieldInfo
from .hilbert import get_intersecting_cpus
from .io_utils import fill_hydro, read_amr
from .particle_handlers import get_particle_handlers


class RAMSESFileSanitizer:
    """A class to handle the different files that can be passed and associated
    safely to a RAMSES output."""

    root_folder = None  # Path | None: path to the root folder
    info_fname = None  # Path | None: path to the info file
    group_name = None  # str | None: name of the first group folder (if any)

    def __init__(self, filename):
        # Make the resolve optional, so that it works with symlinks
        paths_to_try = (Path(filename), Path(filename).resolve())

        self.original_filename = filename

        self.output_dir = None
        self.info_fname = None

        check_functions = (self.test_with_standard_file, self.test_with_folder_name)

        # Loop on both the functions and the tested paths
        for path, check_fun in product(paths_to_try, check_functions):
            ok, output_dir, group_dir, info_fname = check_fun(path)
            if ok:
                break

        # Early exit if the ok flag is False
        if not ok:
            return

        self.root_folder = output_dir
        self.group_name = group_dir.name if group_dir else None
        self.info_fname = info_fname

    def validate(self) -> None:
        # raise a TypeError if self.original_filename is not a valid path
        # we also want to expand '$USER' and '~' because os.path.exists('~') is always False
        filename: str = os.path.expanduser(self.original_filename)

        if not os.path.exists(filename):
            raise FileNotFoundError(rf"No such file or directory '{filename!s}'")
        if self.root_folder is None:
            raise ValueError(
                f"Could not determine output directory from '{filename!s}'\n"
                f"Expected a directory name of form {OUTPUT_DIR_EXP!r} "
                "containing an info_*.txt file and amr_* files."
            )

        # This last case is (erroneously ?) marked as unreachable by mypy
        # If/when this bug is fixed upstream, mypy will warn that the unused
        # 'type: ignore' comment can be removed
        if self.info_fname is None:  # type: ignore [unreachable]
            raise ValueError(f"Failed to detect info file from '{filename!s}'")

    @property
    def is_valid(self) -> bool:
        try:
            self.validate()
        except (TypeError, FileNotFoundError, ValueError):
            return False
        else:
            return True

    @staticmethod
    def check_standard_files(folder, iout):
        """Return True if the folder contains an amr file and the info file."""
        # Check that the "amr_" and "info_" files exist
        ok = (folder / f"amr_{iout}.out00001").is_file()
        ok &= (folder / f"info_{iout}.txt").is_file()
        return ok

    @staticmethod
    def _match_output_and_group(
        path: Path,
    ) -> tuple[Path, Path | None, str | None]:
        # Make sure we work with a directory of the form `output_XXXXX`
        for p in (path, path.parent):
            match = OUTPUT_DIR_RE.match(p.name)

            if match:
                path = p
                break

        if match is None:
            return path, None, None

        iout = match.group(1)

        # See whether a folder named `group_YYYYY` exists
        group_dir = path / "group_00001"
        if group_dir.is_dir():
            return path, group_dir, iout
        else:
            return path, None, iout

    @classmethod
    def test_with_folder_name(
        cls, output_dir: Path
    ) -> tuple[bool, Path | None, Path | None, Path | None]:
        output_dir, group_dir, iout = cls._match_output_and_group(output_dir)
        ok = output_dir.is_dir() and iout is not None

        info_fname: Path | None

        if ok:
            parent_dir = group_dir or output_dir
            ok &= cls.check_standard_files(parent_dir, iout)
            info_fname = parent_dir / f"info_{iout}.txt"
        else:
            info_fname = None

        return ok, output_dir, group_dir, info_fname

    @classmethod
    def test_with_standard_file(
        cls, filename: Path
    ) -> tuple[bool, Path | None, Path | None, Path | None]:
        output_dir, group_dir, iout = cls._match_output_and_group(filename.parent)
        ok = (
            filename.is_file()
            and STANDARD_FILE_RE.match(filename.name) is not None
            and iout is not None
        )

        info_fname: Path | None

        if ok:
            parent_dir = group_dir or output_dir
            ok &= cls.check_standard_files(parent_dir, iout)
            info_fname = parent_dir / f"info_{iout}.txt"
        else:
            info_fname = None

        return ok, output_dir, group_dir, info_fname


class RAMSESDomainFile:
    _last_mask = None
    _last_selector_id = None

    _hydro_offset = None
    _level_count = None

    _oct_handler_initialized = False
    _amr_header_initialized = False

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
        self.max_level = self.ds.parameters["levelmax"] - self.ds.parameters["levelmin"]

        # Autodetect field files
        field_handlers = [FH(self) for FH in get_field_handlers() if FH.any_exist(ds)]
        self.field_handlers = field_handlers
        for fh in field_handlers:
            mylog.debug("Detected fluid type %s in domain_id=%s", fh.ftype, domain_id)
            fh.detect_fields(ds)
            # self._add_ftype(fh.ftype)

        # Autodetect particle files
        particle_handlers = [
            PH(self) for PH in get_particle_handlers() if PH.any_exist(ds)
        ]
        self.particle_handlers = particle_handlers

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
        if hasattr(self, "_amr_file") and not self._amr_file.close:
            self._amr_file.seek(0)
            return self._amr_file

        f = fpu(self.amr_fn)
        self._amr_file = f
        f.seek(0)
        return f

    def _read_amr_header(self):
        if self._amr_header_initialized:
            return
        hvals = {}
        with self.amr_file as f:
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
                self._ngridbound = f.read_vector("i").astype("int64")
            else:
                self._ngridbound = np.zeros(hvals["nlevelmax"], dtype="int64")
            _free_mem = f.read_attrs((("free_mem", 5, "i"),))
            _ordering = f.read_vector("c")
            f.skip(4)
            # Now we're at the tree itself
            # Now we iterate over each level and each CPU.
            position = f.tell()

        self._amr_header = hvals
        self._amr_offset = position

        # The maximum effective level is the deepest level
        # that has a non-zero number of octs
        nocts_to_this_level = hvals["numbl"].sum(axis=1).cumsum()
        self._max_level = (
            np.argwhere(nocts_to_this_level == nocts_to_this_level[-1])[0][0]
            - self.ds.parameters["levelmin"]
            + 1
        )

        # update levelmax
        force_max_level, convention = self.ds._force_max_level
        if convention == "yt":
            force_max_level += self.ds.min_level + 1
        self._amr_header["nlevelmax"] = min(
            force_max_level, self._amr_header["nlevelmax"]
        )
        self._local_oct_count = hvals["numbl"][
            self.ds.min_level :, self.domain_id - 1
        ].sum()
        imin, imax = self.ds.min_level, self._amr_header["nlevelmax"]
        self._total_oct_count = hvals["numbl"][imin:imax, :].sum(axis=0)

        self._amr_header_initialized = True

    @property
    def ngridbound(self):
        self._read_amr_header()
        return self._ngridbound

    @property
    def amr_offset(self):
        self._read_amr_header()
        return self._amr_offset

    @property
    def max_level(self):
        self._read_amr_header()
        return self._max_level

    @max_level.setter
    def max_level(self, value):
        self._max_level = value

    @property
    def total_oct_count(self):
        self._read_amr_header()
        return self._total_oct_count

    @property
    def local_oct_count(self):
        self._read_amr_header()
        return self._local_oct_count

    @property
    def amr_header(self):
        self._read_amr_header()
        return self._amr_header

    @cached_property
    def oct_handler(self):
        """Open the oct file, read in octs level-by-level.
        For each oct, only the position, index, level and domain
        are needed - its position in the octree is found automatically.
        The most important is finding all the information to feed
        oct_handler.add
        """
        self._read_amr_header()
        oct_handler = RAMSESOctreeContainer(
            self.ds.domain_dimensions / 2,
            self.ds.domain_left_edge,
            self.ds.domain_right_edge,
        )
        root_nodes = self.amr_header["numbl"][self.ds.min_level, :].sum()
        oct_handler.allocate_domains(self.total_oct_count, root_nodes)
        mylog.debug(
            "Reading domain AMR % 4i (%0.3e, %0.3e)",
            self.domain_id,
            self.total_oct_count.sum(),
            self.ngridbound.sum(),
        )

        with self.amr_file as f:
            f.seek(self.amr_offset)

            min_level = self.ds.min_level
            max_level = read_amr(
                f, self.amr_header, self.ngridbound, min_level, oct_handler
            )

            oct_handler.finalize()

        new_max_level = max_level
        if new_max_level > self.max_level:
            raise RuntimeError(
                f"The maximum level detected in the AMR file ({new_max_level}) "
                f" does not match the expected number {self.max_level}."
            )
        self.max_level = new_max_level
        self._oct_handler_initialized = True

        return oct_handler

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
        num_zones=2,
        num_ghost_zones=0,
        base_grid=None,
    ):
        super().__init__(base_region, domain, ds, num_zones, num_ghost_zones)

        self._base_grid = base_grid

        if num_ghost_zones > 0:
            if not all(ds.periodicity):
                mylog.warning(
                    "Ghost zones will wrongly assume the domain to be periodic."
                )
            # Create a base domain *with no self._base_domain.fwidth
            base_domain = RAMSESDomainSubset(ds.all_data(), domain, ds, num_zones)
            self._base_domain = base_domain
        elif num_ghost_zones < 0:
            raise RuntimeError(
                "Cannot initialize a domain subset with a negative number "
                f"of ghost zones, was called with {num_ghost_zones=}"
            )

    @property
    def oct_handler(self):
        return self.domain.oct_handler

    def _fill_no_ghostzones(self, fd, fields, selector, file_handler):
        ndim = self.ds.dimensionality
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        oct_handler = self.oct_handler
        all_fields = [f for ft, f in file_handler.field_list]
        fields = [f for ft, f in fields]
        data = {}
        cell_count = selector.count_oct_cells(self.oct_handler, self.domain_id)

        # Initializing data container
        for field in fields:
            data[field] = np.zeros(cell_count, "float64")

        # Do an early exit if the cell count is null
        if cell_count == 0:
            return data

        level_inds, cell_inds, file_inds = self.oct_handler.file_index_octs(
            selector, self.domain_id, cell_count
        )

        cpu_list = [self.domain_id - 1]
        fill_hydro(
            fd,
            file_handler.offset,
            file_handler.level_count,
            cpu_list,
            level_inds,
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
            selector.count_octs(self.oct_handler, self.domain_id) * self.nz**ndim
        )

        # Initializing data container
        for field in fields:
            tr[field] = np.zeros(cell_count, "float64")

        # Do an early exit if the cell count is null
        if cell_count == 0:
            return tr

        gz_cache = getattr(self, "_ghost_zone_cache", None)
        if gz_cache:
            level_inds, cell_inds, file_inds, domain_inds = gz_cache
        else:
            gz_cache = (
                level_inds,
                cell_inds,
                file_inds,
                domain_inds,
            ) = self.oct_handler.file_index_octs_with_ghost_zones(
                selector, self.domain_id, cell_count, self._num_ghost_zones
            )
            self._ghost_zone_cache = gz_cache

        cpu_list = list(range(ncpu))
        fill_hydro(
            fd,
            file_handler.offset,
            file_handler.level_count,
            cpu_list,
            level_inds,
            cell_inds,
            file_inds,
            ndim,
            all_fields,
            fields,
            tr,
            oct_handler,
            domain_inds=domain_inds,
        )
        return tr

    @property
    def fwidth(self):
        fwidth = super().fwidth
        if self._num_ghost_zones > 0:
            fwidth = fwidth.reshape(-1, 8, 3)
            n_oct = fwidth.shape[0]
            # new_fwidth contains the fwidth of the oct+ghost zones
            # this is a constant array in each oct, so we simply copy
            # the oct value using numpy fancy-indexing
            new_fwidth = np.zeros((n_oct, self.nz**3, 3), dtype=fwidth.dtype)
            new_fwidth[:, :, :] = fwidth[:, 0:1, :]
            fwidth = new_fwidth.reshape(-1, 3)
        return fwidth

    @property
    def fcoords(self):
        num_ghost_zones = self._num_ghost_zones
        if num_ghost_zones == 0:
            return super().fcoords

        oh = self.oct_handler

        indices = oh.fill_index(self.selector).reshape(-1, 8)
        oct_inds, cell_inds = oh.fill_octcellindex_neighbours(
            self.selector, self._num_ghost_zones
        )

        N_per_oct = self.nz**3
        oct_inds = oct_inds.reshape(-1, N_per_oct)
        cell_inds = cell_inds.reshape(-1, N_per_oct)

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
                "%s.retrieve_ghost_zones was called with the "
                "`smoothed` argument set to True. This is not supported, "
                "ignoring it.",
                self,
            )
            smoothed = False

        _subset_with_gz = getattr(self, "_subset_with_gz", {})

        try:
            new_subset = _subset_with_gz[ngz]
            mylog.debug(
                "Reusing previous subset with %s ghost zones for domain %s",
                ngz,
                self.domain_id,
            )
        except KeyError:
            new_subset = RAMSESDomainSubset(
                self.base_region,
                self.domain,
                self.ds,
                num_ghost_zones=ngz,
                base_grid=self,
            )
            _subset_with_gz[ngz] = new_subset

        # Cache the fields
        new_subset.get_data(fields)
        self._subset_with_gz = _subset_with_gz

        return new_subset


class RAMSESIndex(OctreeIndex):
    def __init__(self, ds, dataset_type="ramses"):
        self.fluid_field_list = ds._fields_in_file
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)

        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _initialize_oct_handler(self):
        if self.ds._bbox is not None:
            cpu_list = get_intersecting_cpus(self.dataset, self.dataset._bbox)
        else:
            cpu_list = range(self.dataset["ncpu"])

        self.domains = [RAMSESDomainFile(self.dataset, i + 1) for i in cpu_list]

    @cached_property
    def max_level(self):
        force_max_level, convention = self.ds._force_max_level
        if convention == "yt":
            force_max_level += self.ds.min_level + 1
        return min(force_max_level, max(dom.max_level for dom in self.domains))

    @cached_property
    def num_grids(self):
        return sum(
            dom.local_oct_count
            for dom in self.domains  # + dom.ngridbound.sum()
        )

    def _detect_output_fields(self):
        dsl = set()

        # Get the detected particle fields
        for ph in self.domains[0].particle_handlers:
            dsl.update(set(ph.field_offsets.keys()))

        self.particle_field_list = list(dsl)

        # Get the detected fields
        dsl = set()
        for fh in self.domains[0].field_handlers:
            dsl.update(set(fh.field_list))
        self.fluid_field_list = list(dsl)

        self.field_list = self.particle_field_list + self.fluid_field_list

    def _identify_base_chunk(self, dobj):
        use_fast_hilbert = (
            hasattr(dobj, "get_bbox")
            and self.ds.parameters["ordering type"] == "hilbert"
        )
        if getattr(dobj, "_chunk_info", None) is None:
            if use_fast_hilbert:
                idoms = {
                    idom + 1 for idom in get_intersecting_cpus(self.ds, dobj, factor=3)
                }
                # If the oct handler has been initialized, use it
                domains = []
                for dom in self.domains:
                    # Hilbert indexing is conservative, so reject all those that
                    # aren't in the bbox
                    if dom.domain_id not in idoms:
                        continue
                    # If the domain has its oct handler, refine the selection
                    if dom._oct_handler_initialized and not dom.included(dobj.selector):
                        continue
                    mylog.debug("Identified domain %s", dom.domain_id)

                    domains.append(dom)
                if len(domains) >= 1:
                    mylog.info(
                        "Identified % 5d/% 5d intersecting domains (% 5d through hilbert key indexing)",
                        len(domains),
                        len(self.domains),
                        len(idoms),
                    )
            else:
                domains = [dom for dom in self.domains if dom.included(dobj.selector)]
                if len(domains) >= 1:
                    mylog.info("Identified %s intersecting domains", len(domains))
            base_region = getattr(dobj, "base_region", dobj)

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
        levels = sum(dom.level_count for dom in self.domains)
        desc = {"names": ["numcells", "level"], "formats": ["int64"] * 2}
        max_level = self.dataset.min_level + self.dataset.max_level + 2
        self.level_stats = blankRecordArray(desc, max_level)
        self.level_stats["level"] = list(range(max_level))
        self.level_stats["numcells"] = [0 for i in range(max_level)]
        for level in range(self.dataset.min_level + 1):
            self.level_stats[level + 1]["numcells"] = 2 ** (
                level * self.dataset.dimensionality
            )
        for level in range(levels.shape[1]):
            ncell = levels[:, level].sum()
            self.level_stats[level + self.dataset.min_level + 1]["numcells"] = ncell

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

        header = "{:>3}\t{:>14}\t{:>14}".format("level", "# cells", "# cells^3")
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
            print(f"z = {self.dataset.current_redshift:0.8f}")
        except Exception:
            pass
        print(
            "t = {:0.8e} = {:0.8e} = {:0.8e}".format(
                self.ds.current_time.in_units("code_time"),
                self.ds.current_time.in_units("s"),
                self.ds.current_time.in_units("yr"),
            )
        )
        print("\nSmallest Cell:")
        for item in ("Mpc", "pc", "AU", "cm"):
            print(f"\tWidth: {dx.in_units(item):0.3e}")


class RAMSESDataset(Dataset):
    _index_class = RAMSESIndex
    _field_info_class = RAMSESFieldInfo
    gamma = 1.4  # This will get replaced on hydro_fn open

    # RAMSES-specific parameters
    force_cosmological: bool | None
    _force_max_level: tuple[int, str]
    _bbox: list[list[float]] | None
    _self_shielding: bool | None = None

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
        default_species_fields=None,
        self_shielding=None,
        use_conformal_time=None,
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

        self_shielding:
        If set to True, assume gas is self-shielded above 0.01 mp/cm^3.
        This affects the fields related to cooling and the mean molecular weight.
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

        file_handler = RAMSESFileSanitizer(filename)

        # ensure validation happens even if the class is instantiated
        # directly rather than from yt.load
        file_handler.validate()

        # Sanitize the filename
        info_fname = file_handler.info_fname

        if file_handler.group_name is not None:
            self.num_groups = len(
                [_ for _ in file_handler.root_folder.glob("group_?????") if _.is_dir()]
            )
        else:
            self.num_groups = 0
        self.root_folder = file_handler.root_folder

        Dataset.__init__(
            self,
            info_fname,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
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

        if use_conformal_time is not None:
            self.use_conformal_time = use_conformal_time
        elif self.cosmological_simulation:
            if "rt" in self.fluid_types:
                self.use_conformal_time = False
            else:
                self.use_conformal_time = True
        else:
            self.use_conformal_time = False

        self.storage_filename = storage_filename

        self.self_shielding = self_shielding

    @property
    def self_shielding(self) -> bool:
        if self._self_shielding is not None:
            return self._self_shielding

        # Read namelist.txt file (if any)
        has_namelist = self.read_namelist()

        if not has_namelist:
            self._self_shielding = False
            return self._self_shielding

        nml = self.parameters["namelist"]

        # "self_shielding" is stored in physics_params in older versions of the code
        physics_params = nml.get("physics_params", default={})
        # and in "cooling_params" in more recent ones
        cooling_params = nml.get("cooling_params", default={})

        self_shielding = physics_params.get("self_shielding", False)
        self_shielding |= cooling_params.get("self_shielding", False)

        self._self_shielding = self_shielding
        return self_shielding

    @self_shielding.setter
    def self_shielding(self, value):
        self._self_shielding = value

    @staticmethod
    def _sanitize_max_level(max_level, max_level_convention):
        # NOTE: the initialisation of the dataset class sets
        #       self.min_level _and_ requires force_max_level
        #       to be set, so we cannot convert from to yt/ramses
        #       conventions
        if max_level is None and max_level_convention is None:
            return (2**999, "yt")

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
        super().create_field_info(*args, **kwa)
        # Register particle filters
        if ("io", "particle_family") in self.field_list:
            for fname, value in particle_families.items():

                def loc(val):
                    def closure(pfilter, data):
                        filter = data[pfilter.filtered_type, "particle_family"] == val
                        return filter

                    return closure

                add_particle_filter(
                    fname, loc(value), filtered_type="io", requires=["particle_family"]
                )

            for k in particle_families.keys():
                mylog.info("Adding particle_type: %s", k)
                self.add_particle_filter(f"{k}")

    def __str__(self):
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
        mass_unit = density_unit * length_unit**3
        magnetic_unit = np.sqrt(4 * np.pi * mass_unit / (time_unit**2 * length_unit))
        pressure_unit = density_unit * (length_unit / time_unit) ** 2

        setdefaultattr(self, "density_unit", self.quan(density_unit, "g/cm**3"))
        setdefaultattr(self, "magnetic_unit", self.quan(magnetic_unit, "gauss"))
        setdefaultattr(self, "pressure_unit", self.quan(pressure_unit, "dyne/cm**2"))
        setdefaultattr(self, "time_unit", self.quan(time_unit, "s"))
        setdefaultattr(self, "mass_unit", self.quan(mass_unit, "g"))
        setdefaultattr(
            self, "velocity_unit", self.quan(length_unit, "cm") / self.time_unit
        )
        temperature_unit = self.velocity_unit**2 * mp / kb
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
            line = f.readline().strip()

            if line and "=" in line:
                key, val = (_.strip() for _ in line.split("="))
                rheader[key] = cast(val)
                return key
            else:
                return None

        def cast_a_else_b(cast_a, cast_b):
            def caster(val):
                try:
                    return cast_a(val)
                except ValueError:
                    return cast_b(val)

            return caster

        with open(self.parameter_filename) as f:
            # Standard: first six are ncpu, ndim, levelmin, levelmax, ngridmax, nstep_coarse
            for _ in range(6):
                read_rhs(f, int)
            f.readline()
            # Standard: next 11 are boxlen, time, aexp, h0, omega_m, omega_l, omega_k, omega_b, unit_l, unit_d, unit_t
            for _ in range(11):
                key = read_rhs(f, float)

            # Read non standard extra fields until hitting the ordering type
            while key != "ordering type":
                key = read_rhs(f, cast_a_else_b(float, str))

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
                f"The ordering {rheader['ordering type']} "
                "is not compatible with the `bbox` argument."
            )
        self.parameters.update(rheader)
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_dimensions = np.ones(3, dtype="int32") * 2 ** (self.min_level + 1)
        self.domain_right_edge = np.ones(3, dtype="float64")
        # This is likely not true, but it's not clear
        # how to determine the boundary conditions
        self._periodicity = (True, True, True)

        if self.force_cosmological is not None:
            is_cosmological = self.force_cosmological
        else:
            # These conditions seem to always be true for non-cosmological datasets
            is_cosmological = not (
                rheader["time"] >= 0 and rheader["H0"] == 1 and rheader["aexp"] == 1
            )
        if not is_cosmological:
            self.cosmological_simulation = False
            self.current_redshift = 0
            self.hubble_constant = 0
            self.omega_matter = 0
            self.omega_lambda = 0
        else:
            self.cosmological_simulation = True
            self.current_redshift = (1.0 / rheader["aexp"]) - 1.0
            self.omega_lambda = rheader["omega_l"]
            self.omega_matter = rheader["omega_m"]
            self.hubble_constant = rheader["H0"] / 100.0  # This is H100

        force_max_level, convention = self._force_max_level
        if convention == "yt":
            force_max_level += self.min_level + 1
        self.max_level = min(force_max_level, rheader["levelmax"]) - self.min_level - 1

        if not self.cosmological_simulation:
            self.current_time = self.parameters["time"]
        else:
            aexp_grid = np.geomspace(1e-3, 1, 2_000, endpoint=False)
            z_grid = 1 / aexp_grid - 1
            self.tau_frw = tau_frw(self, z_grid)
            self.t_frw = t_frw(self, z_grid)
            self.current_time = t_frw(self, self.current_redshift).to("Gyr")

        if self.num_groups > 0:
            self.group_size = rheader["ncpu"] // self.num_groups

        self.read_namelist()

    def read_namelist(self) -> bool:
        """Read the namelist.txt file in the output folder, if present"""
        namelist_file = os.path.join(self.root_folder, "namelist.txt")
        if not os.path.exists(namelist_file):
            return False

        try:
            with open(namelist_file) as f:
                nml = f90nml.read(f)
        except ImportError as err:
            mylog.warning(
                "`namelist.txt` file found but missing package f90nml to read it:",
                exc_info=err,
            )
            return False
        except (ValueError, StopIteration, AssertionError) as err:
            # Note: f90nml may raise a StopIteration, a ValueError or an AssertionError if
            # the namelist is not valid.
            mylog.warning(
                "Could not parse `namelist.txt` file as it was malformed:",
                exc_info=err,
            )
            return False

        self.parameters["namelist"] = nml
        return True

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if cls._missing_load_requirements():
            return False

        return RAMSESFileSanitizer(filename).is_valid

import os
import re
import weakref
from pathlib import Path

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.frontends.idefix.dmpfile_io import read_header, read_idefix_dmpfile
from yt.frontends.idefix.inifile_io import read_idefix_inifile
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex

from .fields import IdefixFieldInfo


class IdefixGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dims):
        super(IdefixGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims

    def __repr__(self):
        return "IdefixGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class IdefixHierarchy(GridIndex):
    grid = IdefixGrid

    def __init__(self, ds, dataset_type="idefix"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super(IdefixHierarchy, self).__init__(ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [
            (self.dataset_type, f) for f in self.dataset._detected_field_list
        ]

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 1
        self.max_level = 1

    def _populate_grid_objects(self):
        # the minimal form of this method is
        #
        # for g in self.grids:
        #     g._prepare_grid()
        #     g._setup_dx()
        #
        # This must also set:
        #   g.Children <= list of child grids
        #   g.Parent   <= parent grid
        # This is handled by the frontend because often the children must be identified.
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            g = self.grid(i, self, self.grid_levels.flat[i], self.grid_dimensions[i])
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g


class IdefixDataset(Dataset):
    _index_class = IdefixHierarchy
    _field_info_class = IdefixFieldInfo

    def __init__(
        self, dmpfile, inifile=None, dataset_type="idefix", units_override=None,
    ):
        self.fluid_types += ("idefix",)
        self.inifile = inifile
        super(IdefixDataset, self).__init__(
            dmpfile, dataset_type, units_override=units_override
        )
        self.storage_filename = None

        # idefix does not support grid refinement
        self.refine_by = 1

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        # self.mass_unit = self.quan(1.0, "g")
        # self.time_unit = self.quan(1.0, "s")
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    def _parse_parameter_file(self):
        # first pass in the dmpfile: read everything except large arrays
        fprops, fdata = read_idefix_dmpfile(self.parameter_filename, skip_data=True)
        self._detected_field_list = [k for k in fprops if re.match(r"^V[sc]-", k)]

        # parse the version hash
        header = read_header(self.parameter_filename)
        version_exp = r"v\d+\.\d+[-\w]+"  # version number + git hash
        self.parameters["idefix version"] = re.findall(version_exp, header)[0]

        if self.inifile is not None:
            self.parameters.update(read_idefix_inifile(self.inifile))
            grid_ini = self.parameters["Grid"]
            for ax, vals in grid_ini.items():
                if vals[0] > 1:
                    # more than one block is only relevant for mixing grid spacings,
                    # but only "u" is supported
                    raise ValueError(f"Unsupported block structure for {ax}.")
                if vals[3] != "u":
                    raise ValueError(f"Unsupported grid spacing '{vals[3]}'.")

        # parse the grid
        axes = ("x1", "x2", "x3")
        self.domain_dimensions = np.concatenate([fprops[k][-1] for k in axes])
        self.dimensionality = np.count_nonzero(self.domain_dimensions - 1)

        # note that domain edges parsing is already implemented in a mutli-block
        # supporting fashion even though we specifically error out in case there's more
        # than one block.
        self.domain_left_edge = np.array(
            [fdata[f"xl{idir}"][0] for idir in "123"], dtype="float64"
        )
        self.domain_right_edge = np.array(
            [fdata[f"xr{idir}"][-1] for idir in "123"], dtype="float64"
        )

        self.current_time = fdata["time"]

        self.periodicity = tuple(bool(p) for p in fdata["periodicity"])
        enum_geoms = {1: "cartesian", 2: "cylindrial", 3: "polar", 4: "spherical"}
        self.geometry = enum_geoms[fdata["geometry"]]

        # idefix is never cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    @classmethod
    def _is_valid(self, fn, *args, **kwargs):
        ok = bool(re.match(r"^(dump)\.\d{4}(\.dmp)$", Path(fn).name))
        ok &= "idefix" in read_header(fn).lower()
        return ok

import glob
import os
from collections import defaultdict

import numpy as np

from yt.data_objects.static_output import ParticleDataset, ParticleFile
from yt.frontends.gadget.data_structures import _fix_unit_ordering
from yt.funcs import only_on_root, setdefaultattr
from yt.geometry.particle_geometry_handler import ParticleIndex
from yt.utilities.exceptions import YTException
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import OWLSSubfindFieldInfo


class OWLSSubfindParticleIndex(ParticleIndex):
    chunksize = -1

    def __init__(self, ds, dataset_type):
        super().__init__(ds, dataset_type)

    def _calculate_particle_index_starts(self):
        # Halo indices are not saved in the file, so we must count by hand.
        # File 0 has halos 0 to N_0 - 1, file 1 has halos N_0 to N_0 + N_1 - 1, etc.
        particle_count = defaultdict(int)
        offset_count = 0
        for data_file in self.data_files:
            data_file.index_start = {
                ptype: particle_count[ptype] for ptype in data_file.total_particles
            }
            data_file.offset_start = offset_count
            for ptype in data_file.total_particles:
                particle_count[ptype] += data_file.total_particles[ptype]
            offset_count += data_file.total_offset

    def _calculate_file_offset_map(self):
        # After the FOF  is performed, a load-balancing step redistributes halos
        # and then writes more fields.  Here, for each file, we create a list of
        # files which contain the rest of the redistributed particles.
        ifof = np.array(
            [data_file.total_particles["FOF"] for data_file in self.data_files]
        )
        isub = np.array([data_file.total_offset for data_file in self.data_files])
        subend = isub.cumsum()
        fofend = ifof.cumsum()
        istart = np.digitize(fofend - ifof, subend - isub) - 1
        iend = np.clip(np.digitize(fofend, subend), 0, ifof.size - 2)
        for i, data_file in enumerate(self.data_files):
            data_file.offset_files = self.data_files[istart[i] : iend[i] + 1]

    def _detect_output_fields(self):
        # TODO: Add additional fields
        self._calculate_particle_index_starts()
        self._calculate_file_offset_map()
        dsl = []
        units = {}
        for dom in self.data_files:
            fl, _units = self.io._identify_fields(dom)
            units.update(_units)
            for f in fl:
                if f not in dsl:
                    dsl.append(f)
        self.field_list = dsl
        ds = self.dataset
        ds.particle_types = tuple({pt for pt, ds in dsl})
        # This is an attribute that means these particle types *actually*
        # exist.  As in, they are real, in the dataset.
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types


class OWLSSubfindHDF5File(ParticleFile):
    def __init__(self, ds, io, filename, file_id, bounds):
        super().__init__(ds, io, filename, file_id, bounds)
        with h5py.File(filename, mode="r") as f:
            self.header = {field: f.attrs[field] for field in f.attrs.keys()}


class OWLSSubfindDataset(ParticleDataset):
    _index_class = OWLSSubfindParticleIndex
    _file_class = OWLSSubfindHDF5File
    _field_info_class = OWLSSubfindFieldInfo
    _suffix = ".hdf5"

    def __init__(
        self,
        filename,
        dataset_type="subfind_hdf5",
        index_order=None,
        index_filename=None,
        units_override=None,
        unit_system="cgs",
    ):
        super().__init__(
            filename,
            dataset_type,
            index_order=index_order,
            index_filename=index_filename,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _parse_parameter_file(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        hvals = {}
        hvals.update((str(k), v) for k, v in handle["/Header"].attrs.items())
        hvals["NumFiles"] = hvals["NumFilesPerSnapshot"]
        hvals["Massarr"] = hvals["MassTable"]

        self.dimensionality = 3
        self.refine_by = 2

        # Set standard values
        if "Time_GYR" in hvals:
            self.current_time = self.quan(hvals["Time_GYR"], "Gyr")
        self.domain_left_edge = np.zeros(3, "float64")
        self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]
        self.domain_dimensions = np.ones(3, "int32")
        self.cosmological_simulation = 1
        self._periodicity = (True, True, True)
        self.current_redshift = hvals["Redshift"]
        self.omega_lambda = hvals["OmegaLambda"]
        self.omega_matter = hvals["Omega0"]
        self.hubble_constant = hvals["HubbleParam"]
        self.parameters = hvals
        prefix = os.path.abspath(
            os.path.join(
                os.path.dirname(self.parameter_filename),
                os.path.basename(self.parameter_filename).split(".", 1)[0],
            )
        )

        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        self.filename_template = f"{prefix}.%(num)i.{suffix}"
        self.file_count = len(glob.glob(prefix + "*" + self._suffix))
        if self.file_count == 0:
            raise YTException(message="No data files found.", ds=self)
        self.particle_types = ("FOF", "SUBFIND")
        self.particle_types_raw = ("FOF", "SUBFIND")

        # To avoid having to open files twice
        self._unit_base = {}
        self._unit_base.update((str(k), v) for k, v in handle["/Units"].attrs.items())
        handle.close()

    def _set_code_unit_attributes(self):
        # Set a sane default for cosmological simulations.
        if self._unit_base is None and self.cosmological_simulation == 1:
            only_on_root(mylog.info, "Assuming length units are in Mpc/h (comoving)")
            self._unit_base = dict(length=(1.0, "Mpccm/h"))
        # The other same defaults we will use from the standard Gadget
        # defaults.
        unit_base = self._unit_base or {}

        if "length" in unit_base:
            length_unit = unit_base["length"]
        elif "UnitLength_in_cm" in unit_base:
            if self.cosmological_simulation == 0:
                length_unit = (unit_base["UnitLength_in_cm"], "cm")
            else:
                length_unit = (unit_base["UnitLength_in_cm"], "cmcm/h")
        else:
            raise RuntimeError
        length_unit = _fix_unit_ordering(length_unit)
        setdefaultattr(self, "length_unit", self.quan(length_unit[0], length_unit[1]))

        if "velocity" in unit_base:
            velocity_unit = unit_base["velocity"]
        elif "UnitVelocity_in_cm_per_s" in unit_base:
            velocity_unit = (unit_base["UnitVelocity_in_cm_per_s"], "cm/s")
        else:
            velocity_unit = (1e5, "cm/s  * sqrt(a)")
        velocity_unit = _fix_unit_ordering(velocity_unit)
        setdefaultattr(
            self, "velocity_unit", self.quan(velocity_unit[0], velocity_unit[1])
        )

        # We set hubble_constant = 1.0 for non-cosmology, so this is safe.
        # Default to 1e10 Msun/h if mass is not specified.
        if "mass" in unit_base:
            mass_unit = unit_base["mass"]
        elif "UnitMass_in_g" in unit_base:
            if self.cosmological_simulation == 0:
                mass_unit = (unit_base["UnitMass_in_g"], "g")
            else:
                mass_unit = (unit_base["UnitMass_in_g"], "g/h")
        else:
            # Sane default
            mass_unit = (1.0, "1e10*Msun/h")
        mass_unit = _fix_unit_ordering(mass_unit)
        setdefaultattr(self, "mass_unit", self.quan(mass_unit[0], mass_unit[1]))

        if "time" in unit_base:
            time_unit = unit_base["time"]
        elif "UnitTime_in_s" in unit_base:
            time_unit = (unit_base["UnitTime_in_s"], "s")
        else:
            tu = (self.length_unit / self.velocity_unit).to("yr/h")
            time_unit = (tu.d, tu.units)
        setdefaultattr(self, "time_unit", self.quan(time_unit[0], time_unit[1]))

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Constants", "Header", "Parameters", "Units", "FOF"]
        veto_groups = []
        valid = True
        try:
            fh = h5py.File(filename, mode="r")
            valid = all(ng in fh["/"] for ng in need_groups) and not any(
                vg in fh["/"] for vg in veto_groups
            )
            fh.close()
        except Exception:
            valid = False
            pass
        return valid

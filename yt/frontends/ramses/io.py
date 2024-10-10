from collections import defaultdict
from functools import lru_cache
from typing import TYPE_CHECKING, Union

import numpy as np
from unyt import unyt_array

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.frontends.ramses.definitions import VAR_DESC_RE, VERSION_RE
from yt.utilities.cython_fortran_utils import FortranFile
from yt.utilities.exceptions import (
    YTFieldTypeNotFound,
    YTFileNotParseable,
    YTParticleOutputFormatNotImplemented,
)
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

if TYPE_CHECKING:
    import os


def convert_ramses_ages(ds, conformal_ages):
    issue_deprecation_warning(
        msg=(
            "The `convert_ramses_ages' function is deprecated. It should be replaced "
            "by the `convert_ramses_conformal_time_to_physical_age' function."
        ),
        stacklevel=3,
        since="4.0.3",
    )
    return convert_ramses_conformal_time_to_physical_time(ds, conformal_ages)


def convert_ramses_conformal_time_to_physical_time(
    ds, conformal_time: np.ndarray
) -> unyt_array:
    """
    Convert conformal times (as defined in RAMSES) to physical times.

    Arguments
    ---------
    ds : RAMSESDataset
        The RAMSES dataset to use for the conversion
    conformal_time : np.ndarray
        The conformal time as read from disk

    Returns
    -------
    physical_age : np.ndarray
        The physical age in code units
    """
    h0 = ds.hubble_constant
    tau_bins = ds.tau_frw * h0
    t_bins = ds.t_frw

    min_time = 0
    max_time = ds.current_time.to(t_bins.units)

    return ds.arr(
        np.clip(
            np.interp(
                conformal_time,
                tau_bins,
                t_bins.value,
                right=max_time,
                left=min_time,
            ),
            min_time,
            max_time.value,
        ),
        t_bins.units,
    )


def _ramses_particle_binary_file_handler(particle_handler, subset, fields, count):
    """General file handler for binary file, called by _read_particle_subset

    Parameters
    ----------
    particle : ``ParticleFileHandler``
        the particle class we want to read
    subset: ``RAMSESDomainSubset``
        A RAMSES domain subset object
    fields: list of tuple
        The fields to read
    count: integer
        The number of elements to count
    """
    tr = {}
    ds = subset.domain.ds
    foffsets = particle_handler.field_offsets
    fname = particle_handler.fname
    data_types = particle_handler.field_types
    with FortranFile(fname) as fd:
        # We do *all* conversion into boxlen here.
        # This means that no other conversions need to be applied to convert
        # positions into the same domain as the octs themselves.
        for field in sorted(fields, key=lambda a: foffsets[a]):
            if count == 0:
                tr[field] = np.empty(0, dtype=data_types[field])
                continue
            # Sentinel value: -1 means we don't have this field
            if foffsets[field] == -1:
                tr[field] = np.empty(count, dtype=data_types[field])
            else:
                fd.seek(foffsets[field])
                dt = data_types[field]
                tr[field] = fd.read_vector(dt)
            if field[1].startswith("particle_position"):
                np.divide(tr[field], ds["boxlen"], tr[field])

            # Hand over to field handler for special cases, like particle_birth_times
            particle_handler.handle_field(field, tr)

    return tr


def _ramses_particle_csv_file_handler(particle_handler, subset, fields, count):
    """General file handler for csv file, called by _read_particle_subset

    Parameters
    ----------
    particle: ``ParticleFileHandler``
        the particle class we want to read
    subset: ``RAMSESDomainSubset``
        A RAMSES domain subset object
    fields: list of tuple
        The fields to read
    count: integer
        The number of elements to count
    """
    from yt.utilities.on_demand_imports import _pandas as pd

    tr = {}
    ds = subset.domain.ds
    foffsets = particle_handler.field_offsets
    fname = particle_handler.fname

    list_field_ind = [
        (field, foffsets[field]) for field in sorted(fields, key=lambda a: foffsets[a])
    ]

    # read only selected fields
    dat = pd.read_csv(
        fname,
        delimiter=",",
        usecols=[ind for _field, ind in list_field_ind],
        skiprows=2,
        header=None,
    )

    for field, ind in list_field_ind:
        tr[field] = dat[ind].to_numpy()
        if field[1].startswith("particle_position"):
            np.divide(tr[field], ds["boxlen"], tr[field])

        particle_handler.handle_field(field, tr)

    return tr


class IOHandlerRAMSES(BaseIOHandler):
    _dataset_type = "ramses"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        tr = defaultdict(list)

        # Set of field types
        ftypes = {f[0] for f in fields}
        for chunk in chunks:
            # Gather fields by type to minimize i/o operations
            for ft in ftypes:
                # Get all the fields of the same type
                field_subs = list(filter(lambda f, ft=ft: f[0] == ft, fields))

                # Loop over subsets
                for subset in chunk.objs:
                    fname = None
                    for fh in subset.domain.field_handlers:
                        if fh.ftype == ft:
                            file_handler = fh
                            fname = fh.fname
                            break

                    if fname is None:
                        raise YTFieldTypeNotFound(ft)

                    # Now we read the entire thing
                    with FortranFile(fname) as fd:
                        # This contains the boundary information, so we skim through
                        # and pick off the right vectors
                        rv = subset.fill(fd, field_subs, selector, file_handler)
                    for ft, f in field_subs:
                        d = rv.pop(f)
                        if d.size == 0:
                            continue
                        mylog.debug(
                            "Filling %s with %s (%0.3e %0.3e) (%s zones)",
                            f,
                            d.size,
                            d.min(),
                            d.max(),
                            d.size,
                        )
                        tr[ft, f].append(d)
        d = {}
        for field in fields:
            tmp = tr.pop(field, None)
            d[field] = np.concatenate(tmp) if tmp else np.empty(0, dtype="d")

        return d

    def _read_particle_coords(self, chunks, ptf):
        pn = "particle_position_%s"
        fields = [
            (ptype, f"particle_position_{ax}")
            for ptype, field_list in ptf.items()
            for ax in "xyz"
        ]
        for chunk in chunks:
            for subset in chunk.objs:
                rv = self._read_particle_subset(subset, fields)
                for ptype in sorted(ptf):
                    yield (
                        ptype,
                        (
                            rv[ptype, pn % "x"],
                            rv[ptype, pn % "y"],
                            rv[ptype, pn % "z"],
                        ),
                        0.0,
                    )

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = "particle_position_%s"
        chunks = list(chunks)
        fields = [
            (ptype, fname) for ptype, field_list in ptf.items() for fname in field_list
        ]

        for ptype, field_list in sorted(ptf.items()):
            for ax in "xyz":
                if pn % ax not in field_list:
                    fields.append((ptype, pn % ax))

            if ptype == "sink_csv":
                subset = chunks[0].objs[0]
                rv = self._read_particle_subset(subset, fields)
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (np.asarray(rv[ptype, pn % ax], "=f8") for ax in "xyz")
                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        mask = []
                    for field in field_list:
                        data = np.asarray(rv.pop((ptype, field))[mask], "=f8")
                        yield (ptype, field), data

            else:
                for chunk in chunks:
                    for subset in chunk.objs:
                        rv = self._read_particle_subset(subset, fields)
                        for ptype, field_list in sorted(ptf.items()):
                            x, y, z = (
                                np.asarray(rv[ptype, pn % ax], "=f8") for ax in "xyz"
                            )
                            mask = selector.select_points(x, y, z, 0.0)
                            if mask is None:
                                mask = []
                            for field in field_list:
                                data = np.asarray(rv.pop((ptype, field))[mask], "=f8")
                                yield (ptype, field), data

    def _read_particle_subset(self, subset, fields):
        """Read the particle files."""
        tr = {}

        # Sequential read depending on particle type
        for ptype in {f[0] for f in fields}:
            # Select relevant files
            subs_fields = filter(lambda f, ptype=ptype: f[0] == ptype, fields)

            ok = False
            for ph in subset.domain.particle_handlers:
                if ph.ptype == ptype:
                    ok = True
                    count = ph.local_particle_count
                    break
            if not ok:
                raise YTFieldTypeNotFound(ptype)

            tr.update(ph.reader(subset, subs_fields, count))

        return tr


@lru_cache
def _read_part_binary_file_descriptor(fname: Union[str, "os.PathLike[str]"]):
    """
    Read a file descriptor and returns the array of the fields found.
    """

    # Mapping
    mapping_list = [
        ("position_x", "particle_position_x"),
        ("position_y", "particle_position_y"),
        ("position_z", "particle_position_z"),
        ("velocity_x", "particle_velocity_x"),
        ("velocity_y", "particle_velocity_y"),
        ("velocity_z", "particle_velocity_z"),
        ("mass", "particle_mass"),
        ("identity", "particle_identity"),
        ("levelp", "particle_level"),
        ("family", "particle_family"),
        ("tag", "particle_tag"),
    ]
    # Convert to dictionary
    mapping = dict(mapping_list)

    with open(fname) as f:
        line = f.readline()
        tmp = VERSION_RE.match(line)
        mylog.debug("Reading part file descriptor %s.", fname)
        if not tmp:
            raise YTParticleOutputFormatNotImplemented()

        version = int(tmp.group(1))

        if version == 1:
            # Skip one line (containing the headers)
            line = f.readline()
            fields = []
            for i, line in enumerate(f.readlines()):
                tmp = VAR_DESC_RE.match(line)
                if not tmp:
                    raise YTFileNotParseable(fname, i + 1)

                # ivar = tmp.group(1)
                varname = tmp.group(2)
                dtype = tmp.group(3)

                if varname in mapping:
                    varname = mapping[varname]
                else:
                    varname = f"particle_{varname}"

                fields.append((varname, dtype))
        else:
            raise YTParticleOutputFormatNotImplemented()

    return fields


@lru_cache
def _read_part_csv_file_descriptor(fname: Union[str, "os.PathLike[str]"]):
    """
    Read the file from the csv sink particles output.
    """
    from yt.utilities.on_demand_imports import _pandas as pd

    # Fields name from the default csv RAMSES sink algorithm in the yt default convention
    mapping = {
        " # id": "particle_identifier",
        "msink": "particle_mass",
        "x": "particle_position_x",
        "y": "particle_position_y",
        "z": "particle_position_z",
        "vx": "particle_velocity_x",
        "vy": "particle_velocity_y",
        "vz": "particle_velocity_z",
        "lx": "particle_angular_momentum_x",
        "ly": "particle_angular_momentum_y",
        "lz": "particle_angular_momentum_z",
        "tform": "particle_formation_time",
        "acc_rate": "particle_accretion_rate",
        "del_mass": "particle_delta_mass",
        "rho_gas": "particle_rho_gas",
        "cs**2": "particle_sound_speed",
        "etherm": "particle_etherm",
        "vx_gas": "particle_velocity_x_gas",
        "vy_gas": "particle_velocity_y_gas",
        "vz_gas": "particle_velocity_z_gas",
        "mbh": "particle_mass_bh",
        "level": "particle_level",
        "rsink_star": "particle_radius_star",
    }

    # read the all file to get the number of particle
    dat = pd.read_csv(fname, delimiter=",")
    fields = []
    local_particle_count = len(dat)

    for varname in dat.columns:
        if varname in mapping:
            varname = mapping[varname]
        else:
            varname = f"particle_{varname}"

        fields.append(varname)

    return fields, local_particle_count


@lru_cache
def _read_fluid_file_descriptor(fname: Union[str, "os.PathLike[str]"], *, prefix: str):
    """
    Read a file descriptor and returns the array of the fields found.
    """

    # Mapping
    mapping_list = [
        ("density", "Density"),
        ("velocity_x", "x-velocity"),
        ("velocity_y", "y-velocity"),
        ("velocity_z", "z-velocity"),
        ("pressure", "Pressure"),
        ("metallicity", "Metallicity"),
        # Add mapping for ionized species
        # Note: we expect internally that these names use the HII, HeII,
        #       HeIII, ... convention for historical reasons. So we need to map
        #       the names read from `hydro_file_descriptor.txt` to this
        #       convention.
        # This will create fields like ("ramses", "HII") which are mapped
        # to ("gas", "H_p1_fraction") in fields.py
        ("H_p1_fraction", "HII"),
        ("He_p1_fraction", "HeII"),
        ("He_p2_fraction", "HeIII"),
        # Photon fluxes / densities are stored as `photon_density_XX`, so
        # only 100 photon bands can be stored with this format. Let's be
        # conservative and support up to 100 bands.
        *[(f"photon_density_{i:02d}", f"Photon_density_{i:d}") for i in range(100)],
        *[
            (f"photon_flux_{i:02d}_{dim}", f"Photon_flux_{dim}_{i:d}")
            for i in range(100)
            for dim in "xyz"
        ],
    ]

    # Add mapping for magnetic fields
    mapping_list += [
        (key, key)
        for key in (
            f"B_{dim}_{side}" for side in ["left", "right"] for dim in ["x", "y", "z"]
        )
    ]

    # Convert to dictionary
    mapping = dict(mapping_list)

    with open(fname) as f:
        line = f.readline()
        tmp = VERSION_RE.match(line)
        mylog.debug("Reading fluid file descriptor %s.", fname)
        if not tmp:
            return []

        version = int(tmp.group(1))

        if version == 1:
            # Skip one line (containing the headers)
            line = f.readline()
            fields = []
            for i, line in enumerate(f.readlines()):
                tmp = VAR_DESC_RE.match(line)
                if not tmp:
                    raise YTFileNotParseable(fname, i + 1)

                # ivar = tmp.group(1)
                varname = tmp.group(2)
                dtype = tmp.group(3)

                if varname in mapping:
                    varname = mapping[varname]
                else:
                    varname = f"{prefix}_{varname}"

                fields.append((varname, dtype))
        else:
            mylog.error("Version %s", version)
            raise YTParticleOutputFormatNotImplemented()

    return fields

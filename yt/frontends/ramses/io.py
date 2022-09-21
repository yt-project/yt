import os
from collections import defaultdict
from typing import Union

import numpy as np

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
from yt.utilities.physical_ratios import cm_per_km, cm_per_mpc


def convert_ramses_ages(ds, conformal_ages):
    issue_deprecation_warning(
        msg=(
            "The `convert_ramses_ages' function is deprecated. It should be replaced "
            "by the `convert_ramses_conformal_time_to_physical_age' function."
        ),
        since="4.0.2",
    )
    return convert_ramses_conformal_time_to_physical_age(ds, conformal_ages)


def convert_ramses_conformal_time_to_physical_age(
    ds, conformal_time: np.ndarray
) -> np.ndarray:
    """
    Convert conformal times (as defined in RAMSES) to physical age.

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
    tf = ds.t_frw
    dtau = ds.dtau
    tauf = ds.tau_frw
    h100 = ds.hubble_constant
    nOver2 = ds.n_frw / 2
    unit_t = ds.parameters["unit_t"]
    tsim = ds.time_simu

    t_scale = 1.0 / (h100 * 100 * cm_per_km / cm_per_mpc) / unit_t

    # calculate index into lookup table (n_frw elements in
    # lookup table)
    dage = 1 + (10 * conformal_time / dtau)
    dage = np.minimum(dage, nOver2 + (dage - nOver2) / 10.0)
    iage = np.array(dage, dtype=np.int32)

    # linearly interpolate physical times from tf and tauf lookup
    # tables.
    t = tf[iage] * (conformal_time - tauf[iage - 1]) / (tauf[iage] - tauf[iage - 1])
    t = t + (
        tf[iage - 1] * (conformal_time - tauf[iage]) / (tauf[iage - 1] - tauf[iage])
    )
    return (tsim - t) * t_scale


def _ramses_particle_file_handler(fname, foffsets, data_types, subset, fields, count):
    """General file handler, called by _read_particle_subset

    Parameters
    ----------
    fname : string
        filename to read from
    foffsets: dict
        Offsets in file of the fields
    data_types: dict
        Data type of the fields
    subset: ``RAMSESDomainSubset``
        A RAMSES domain subset object
    fields: list of tuple
        The fields to read
    count: integer
        The number of elements to count
    """
    tr = {}
    ds = subset.domain.ds
    current_time = ds.current_time.in_units("code_time").v
    with FortranFile(fname) as fd:
        # We do *all* conversion into boxlen here.
        # This means that no other conversions need to be applied to convert
        # positions into the same domain as the octs themselves.
        for field in sorted(fields, key=lambda a: foffsets[a]):
            if count == 0:
                tr[field] = np.empty(0, dtype=data_types[field])
                continue
            fd.seek(foffsets[field])
            dt = data_types[field]
            tr[field] = fd.read_vector(dt)
            if field[1].startswith("particle_position"):
                np.divide(tr[field], ds["boxlen"], tr[field])
            if ds.cosmological_simulation and field[1] == "particle_birth_time":
                conformal_time = tr[field]
                physical_age = convert_ramses_conformal_time_to_physical_age(
                    ds, conformal_time
                )
                tr[field] = current_time - physical_age
                # arbitrarily set particles with zero conformal_age to zero
                # particle_age. This corresponds to DM particles.
                tr[field][conformal_time == 0] = 0
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
                        mylog.debug(
                            "Filling %s with %s (%0.3e %0.3e) (%s zones)",
                            f,
                            d.size,
                            d.min(),
                            d.max(),
                            d.size,
                        )
                        tr[(ft, f)].append(d)
        d = {}
        for field in fields:
            d[field] = np.concatenate(tr.pop(field))

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
                    yield ptype, (
                        rv[ptype, pn % "x"],
                        rv[ptype, pn % "y"],
                        rv[ptype, pn % "z"],
                    ), 0.0

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
        for chunk in chunks:
            for subset in chunk.objs:
                rv = self._read_particle_subset(subset, fields)
                for ptype, field_list in sorted(ptf.items()):
                    x, y, z = (np.asarray(rv[ptype, pn % ax], "=f8") for ax in "xyz")
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
                    fname = ph.fname
                    foffsets = ph.field_offsets
                    data_types = ph.field_types
                    ok = True
                    count = ph.local_particle_count
                    break
            if not ok:
                raise YTFieldTypeNotFound(ptype)

            cosmo = self.ds.cosmological_simulation
            if (ptype, "particle_birth_time") in foffsets and cosmo:
                foffsets[ptype, "conformal_birth_time"] = foffsets[
                    ptype, "particle_birth_time"
                ]
                data_types[ptype, "conformal_birth_time"] = data_types[
                    ptype, "particle_birth_time"
                ]

            tr.update(
                _ramses_particle_file_handler(
                    fname, foffsets, data_types, subset, subs_fields, count=count
                )
            )

        return tr


def _read_part_file_descriptor(fname: Union[str, "os.PathLike[str]"]):
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
    mapping = {k: v for k, v in mapping_list}

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


def _read_fluid_file_descriptor(fname: Union[str, "os.PathLike[str]"]):
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
    ]

    # Add mapping for magnetic fields
    mapping_list += [
        (key, key)
        for key in (
            f"B_{dim}_{side}" for side in ["left", "right"] for dim in ["x", "y", "z"]
        )
    ]

    # Convert to dictionary
    mapping = {k: v for k, v in mapping_list}

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
                    varname = f"hydro_{varname}"

                fields.append((varname, dtype))
        else:
            mylog.error("Version %s", version)
            raise YTParticleOutputFormatNotImplemented()

    return fields

"""
This module gathers all user-facing functions with a `load_` prefix.

"""
import os
import sys
import tarfile
from pathlib import Path
from typing import List, Optional, Tuple
from urllib.parse import urlsplit

import numpy as np
from more_itertools import always_iterable

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.funcs import levenshtein_distance
from yt.sample_data.api import lookup_on_disk_data
from yt.utilities.decompose import decompose_array, get_psize
from yt.utilities.exceptions import (
    YTAmbiguousDataType,
    YTIllDefinedAMR,
    YTSimulationNotIdentified,
    YTUnidentifiedDataType,
)
from yt.utilities.hierarchy_inspection import find_lowest_subclasses
from yt.utilities.lib.misc_utilities import get_box_grids_level
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.object_registries import (
    output_type_registry,
    simulation_time_series_registry,
)
from yt.utilities.on_demand_imports import _pooch as pooch

# --- Loaders for known data formats ---


def load(fn, *args, **kwargs):
    """
    Load a Dataset or DatasetSeries object.
    The data format is automatically discovered, and the exact return type is the
    corresponding subclass of :class:`yt.data_objects.static_output.Dataset`.
    A :class:`yt.data_objects.time_series.DatasetSeries` is created if the first
    argument is a pattern.

    Parameters
    ----------
    fn : str, os.Pathlike, or byte (types supported by os.path.expandusers)
        A path to the data location. This can be a file name, directory name, a glob
        pattern, or a url (for data types that support it).

    Additional arguments, if any, are passed down to the return class.

    Returns
    -------
    :class:`yt.data_objects.static_output.Dataset` object
        If fn is a single path, create a Dataset from the appropriate subclass.

    :class:`yt.data_objects.time_series.DatasetSeries`
        If fn is a glob pattern (i.e. containing wildcards '[]?!*'), create a series.

    Raises
    ------
    FileNotFoundError
        If fn does not match any existing file or directory.

    yt.utilities.exceptions.YTUnidentifiedDataType
        If fn matches existing files or directories with undetermined format.

    yt.utilities.exceptions.YTAmbiguousDataType
        If the data format matches more than one class of similar specilization levels.
    """
    fn = os.path.expanduser(fn)

    if any(wildcard in fn for wildcard in "[]?!*"):
        from yt.data_objects.time_series import DatasetSeries

        return DatasetSeries(fn, *args, **kwargs)

    # This will raise FileNotFoundError if the path isn't matched
    # either in the current dir or yt.config.ytcfg['data_dir_directory']
    if not fn.startswith("http"):
        fn = str(lookup_on_disk_data(fn))

    candidates = []
    for cls in output_type_registry.values():
        if cls._is_valid(fn, *args, **kwargs):
            candidates.append(cls)

    # Find only the lowest subclasses, i.e. most specialised front ends
    candidates = find_lowest_subclasses(candidates)

    if len(candidates) == 1:
        return candidates[0](fn, *args, **kwargs)

    if len(candidates) > 1:
        raise YTAmbiguousDataType(fn, candidates)

    raise YTUnidentifiedDataType(fn, *args, **kwargs)


def load_simulation(fn, simulation_type, find_outputs=False):
    """
    Load a simulation time series object of the specified simulation type.

    Parameters
    ----------
    fn : str, os.Pathlike, or byte (types supported by os.path.expandusers)
        Name of the data file or directory.

    simulation_type : str
        E.g. 'Enzo'

    find_outputs : bool
        Defaults to False

    Raises
    ------
    FileNotFoundError
        If fn is not found.

    yt.utilities.exceptions.YTSimulationNotIdentified
        If simulation_type is unknown.
    """

    fn = str(lookup_on_disk_data(fn))

    try:
        cls = simulation_time_series_registry[simulation_type]
    except KeyError as e:
        raise YTSimulationNotIdentified(simulation_type) from e

    return cls(fn, find_outputs=find_outputs)


def simulation(fn, simulation_type, find_outputs=False):
    issue_deprecation_warning(
        "yt.simulation is a deprecated alias for yt.load_simulation"
        "and will be removed in a future version of yt.",
        since="4.0.0",
        removal="4.1.0",
    )
    return load_simulation(
        fn=fn, simulation_type=simulation_type, find_outputs=find_outputs
    )


# --- Loaders for generic ("stream") data ---


def load_uniform_grid(
    data,
    domain_dimensions,
    length_unit=None,
    bbox=None,
    nprocs=1,
    sim_time=0.0,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    geometry="cartesian",
    unit_system="cgs",
    default_species_fields=None,
):
    r"""Load a uniform grid of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.

    This should allow a uniform grid of data to be loaded directly into yt and
    analyzed as would any others.  This comes with several caveats:

    * Units will be incorrect unless the unit system is explicitly
      specified.
    * Some functions may behave oddly, and parallelism will be
      disappointing or non-existent in most cases.
    * Particles may be difficult to integrate.

    Particle fields are detected as one-dimensional fields.

    Parameters
    ----------
    data : dict
        This is a dict of numpy arrays or (numpy array, unit spec) tuples.
        The keys are the field names.
    domain_dimensions : array_like
        This is the domain dimensions of the grid
    length_unit : string
        Unit to use for lengths.  Defaults to unitless.
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units specified by length_unit.
        Defaults to a cubic unit-length domain.
    nprocs: integer, optional
        If greater than 1, will create this number of subarrays out of data
    sim_time : float, optional
        The simulation time in seconds
    mass_unit : string
        Unit to use for masses.  Defaults to unitless.
    time_unit : string
        Unit to use for times.  Defaults to unitless.
    velocity_unit : string
        Unit to use for velocities.  Defaults to unitless.
    magnetic_unit : string
        Unit to use for magnetic fields. Defaults to unitless.
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    geometry : string or tuple
        "cartesian", "cylindrical", "polar", "spherical", "geographic" or
        "spectral_cube".  Optionally, a tuple can be provided to specify the
        axis ordering -- for instance, to specify that the axis ordering should
        be z, x, y, this would be: ("cartesian", ("z", "x", "y")).  The same
        can be done for other coordinates, for instance:
        ("spherical", ("theta", "phi", "r")).
    default_species_fields : string, optional
        If set, default species fields are created for H and He which also
        determine the mean molecular weight. Options are "ionized" and "neutral".

    Examples
    --------
    >>> np.random.seed(int(0x4D3D3D3))
    >>> bbox = np.array([[0.0, 1.0], [-1.5, 1.5], [1.0, 2.5]])
    >>> arr = np.random.random((128, 128, 128))
    >>> data = dict(density=arr)
    >>> ds = load_uniform_grid(data, arr.shape, length_unit="cm", bbox=bbox, nprocs=12)
    >>> dd = ds.all_data()
    >>> dd[("gas", "density")]
    unyt_array([0.76017901, 0.96855994, 0.49205428, ..., 0.78798258,
                0.97569432, 0.99453904], 'g/cm**3')
    """
    from yt.frontends.stream.data_structures import (
        StreamDataset,
        StreamDictFieldHandler,
        StreamHandler,
    )
    from yt.frontends.stream.definitions import (
        assign_particle_data,
        process_data,
        set_particle_types,
    )

    domain_dimensions = np.array(domain_dimensions)
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], "float64")
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros(nprocs, dtype="int32").reshape((nprocs, 1))
    # If someone included this throw it away--old API
    if "number_of_particles" in data:
        issue_deprecation_warning(
            "It is no longer necessary to include "
            "the number of particles in the data "
            "dict. The number of particles is "
            "determined from the sizes of the "
            "particle fields.",
            since="4.0.0",
            removal="4.1.0",
        )
        data.pop("number_of_particles")
    # First we fix our field names, apply units to data
    # and check for consistency of field shapes
    field_units, data, number_of_particles = process_data(
        data, grid_dims=tuple(domain_dimensions)
    )

    sfh = StreamDictFieldHandler()

    if number_of_particles > 0:
        particle_types = set_particle_types(data)
        # Used much further below.
        pdata = {"number_of_particles": number_of_particles}
        for key in list(data.keys()):
            if len(data[key].shape) == 1 or key[0] == "io":
                if not isinstance(key, tuple):
                    field = ("io", key)
                    mylog.debug("Reassigning '%s' to '%s'", key, field)
                else:
                    field = key
                sfh._additional_fields += (field,)
                pdata[field] = data.pop(key)
    else:
        particle_types = {}

    if nprocs > 1:
        temp = {}
        new_data = {}
        for key in data.keys():
            psize = get_psize(np.array(data[key].shape), nprocs)
            grid_left_edges, grid_right_edges, shapes, slices = decompose_array(
                data[key].shape, psize, bbox
            )
            grid_dimensions = np.array([shape for shape in shapes], dtype="int32")
            temp[key] = [data[key][slice] for slice in slices]
        for gid in range(nprocs):
            new_data[gid] = {}
            for key in temp.keys():
                new_data[gid].update({key: temp[key][gid]})
        sfh.update(new_data)
        del new_data, temp
    else:
        sfh.update({0: data})
        grid_left_edges = domain_left_edge
        grid_right_edges = domain_right_edge
        grid_dimensions = domain_dimensions.reshape(nprocs, 3).astype("int32")

    if length_unit is None:
        length_unit = "code_length"
    if mass_unit is None:
        mass_unit = "code_mass"
    if time_unit is None:
        time_unit = "code_time"
    if velocity_unit is None:
        velocity_unit = "code_velocity"
    if magnetic_unit is None:
        magnetic_unit = "code_magnetic"

    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype="int64"),
        np.zeros(nprocs, dtype="int64").reshape(nprocs, 1),  # particle count
        np.zeros(nprocs).reshape((nprocs, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "UniformGridData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    if np.all(domain_dimensions[1:] == 1):
        dimensionality = 1
    elif domain_dimensions[2] == 1:
        dimensionality = 2
    else:
        dimensionality = 3
    handler.dimensionality = dimensionality
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamDataset(
        handler,
        geometry=geometry,
        unit_system=unit_system,
        default_species_fields=default_species_fields,
    )

    # Now figure out where the particles go
    if number_of_particles > 0:
        # This will update the stream handler too
        assign_particle_data(sds, pdata, bbox)

    return sds


def load_amr_grids(
    grid_data,
    domain_dimensions,
    bbox=None,
    sim_time=0.0,
    length_unit=None,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    geometry="cartesian",
    refine_by=2,
    unit_system="cgs",
    default_species_fields=None,
):
    r"""Load a set of grids of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.
    This should allow a sequence of grids of varying resolution of data to be
    loaded directly into yt and analyzed as would any others.  This comes with
    several caveats:

    * Units will be incorrect unless the unit system is explicitly specified.
    * Some functions may behave oddly, and parallelism will be
      disappointing or non-existent in most cases.
    * Particles may be difficult to integrate.
    * No consistency checks are performed on the index

    Parameters
    ----------

    grid_data : list of dicts
        This is a list of dicts. Each dict must have entries "left_edge",
        "right_edge", "dimensions", "level", and then any remaining entries are
        assumed to be fields. Field entries must map to an NDArray. The grid_data
        may also include a particle count. If no particle count is supplied, the
        dataset is understood to contain no particles. The grid_data will be
        modified in place and can't be assumed to be static.
    domain_dimensions : array_like
        This is the domain dimensions of the grid
    length_unit : string or float
        Unit to use for lengths.  Defaults to unitless.  If set to be a string, the bbox
        dimensions are assumed to be in the corresponding units.  If set to a float, the
        value is a assumed to be the conversion from bbox dimensions to centimeters.
    mass_unit : string or float
        Unit to use for masses.  Defaults to unitless.
    time_unit : string or float
        Unit to use for times.  Defaults to unitless.
    velocity_unit : string or float
        Unit to use for velocities.  Defaults to unitless.
    magnetic_unit : string or float
        Unit to use for magnetic fields.  Defaults to unitless.
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units specified by length_unit.
        Defaults to a cubic unit-length domain.
    sim_time : float, optional
        The simulation time in seconds
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    geometry : string or tuple
        "cartesian", "cylindrical", "polar", "spherical", "geographic" or
        "spectral_cube".  Optionally, a tuple can be provided to specify the
        axis ordering -- for instance, to specify that the axis ordering should
        be z, x, y, this would be: ("cartesian", ("z", "x", "y")).  The same
        can be done for other coordinates, for instance:
        ("spherical", ("theta", "phi", "r")).
    refine_by : integer or list/array of integers.
        Specifies the refinement ratio between levels.  Defaults to 2.  This
        can be an array, in which case it specifies for each dimension.  For
        instance, this can be used to say that some datasets have refinement of
        1 in one dimension, indicating that they span the full range in that
        dimension.
    default_species_fields : string, optional
        If set, default species fields are created for H and He which also
        determine the mean molecular weight. Options are "ionized" and "neutral".

    Examples
    --------

    >>> grid_data = [
    ...     dict(
    ...         left_edge=[0.0, 0.0, 0.0],
    ...         right_edge=[1.0, 1.0, 1.0],
    ...         level=0,
    ...         dimensions=[32, 32, 32],
    ...         number_of_particles=0,
    ...     ),
    ...     dict(
    ...         left_edge=[0.25, 0.25, 0.25],
    ...         right_edge=[0.75, 0.75, 0.75],
    ...         level=1,
    ...         dimensions=[32, 32, 32],
    ...         number_of_particles=0,
    ...     ),
    ... ]
    ...
    >>> for g in grid_data:
    ...     g[("gas", "density")] = (
    ...         np.random.random(g["dimensions"]) * 2 ** g["level"],
    ...         "g/cm**3",
    ...     )
    ...
    >>> ds = load_amr_grids(grid_data, [32, 32, 32], length_unit=1.0)
    """
    from yt.frontends.stream.data_structures import (
        StreamDataset,
        StreamDictFieldHandler,
        StreamHandler,
    )
    from yt.frontends.stream.definitions import process_data, set_particle_types

    domain_dimensions = np.array(domain_dimensions)
    ngrids = len(grid_data)
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], "float64")
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros((ngrids, 1), dtype="int32")
    grid_left_edges = np.zeros((ngrids, 3), dtype="float64")
    grid_right_edges = np.zeros((ngrids, 3), dtype="float64")
    grid_dimensions = np.zeros((ngrids, 3), dtype="int32")
    number_of_particles = np.zeros((ngrids, 1), dtype="int64")
    parent_ids = np.zeros(ngrids, dtype="int64") - 1
    sfh = StreamDictFieldHandler()
    for i, g in enumerate(grid_data):
        grid_left_edges[i, :] = g.pop("left_edge")
        grid_right_edges[i, :] = g.pop("right_edge")
        grid_dimensions[i, :] = g.pop("dimensions")
        grid_levels[i, :] = g.pop("level")
        # If someone included this throw it away--old API
        if "number_of_particles" in g:
            issue_deprecation_warning(
                "It is no longer necessary to include "
                "the number of particles in the data "
                "dict. The number of particles is "
                "determined from the sizes of the "
                "particle fields.",
                since="4.0.0",
                removal="4.1.0",
            )
            g.pop("number_of_particles")
        field_units, data, n_particles = process_data(
            g, grid_dims=tuple(grid_dimensions[i, :])
        )
        number_of_particles[i, :] = n_particles
        sfh[i] = data

    # We now reconstruct our parent ids, so that our particle assignment can
    # proceed.
    mask = np.empty(ngrids, dtype="int32")
    for gi in range(ngrids):
        get_box_grids_level(
            grid_left_edges[gi, :],
            grid_right_edges[gi, :],
            grid_levels[gi] + 1,
            grid_left_edges,
            grid_right_edges,
            grid_levels,
            mask,
        )
        ids = np.where(mask.astype("bool"))
        for ci in ids:
            parent_ids[ci] = gi

    # Check if the grid structure is properly aligned (bug #1295)
    for lvl in range(grid_levels.min() + 1, grid_levels.max() + 1):
        idx = grid_levels.flatten() == lvl
        dims = domain_dimensions * refine_by ** (lvl - 1)
        for iax, ax in enumerate("xyz"):
            cell_edges = np.linspace(
                domain_left_edge[iax], domain_right_edge[iax], dims[iax], endpoint=False
            )
            if set(grid_left_edges[idx, iax]) - set(cell_edges):
                raise YTIllDefinedAMR(lvl, ax)

    if length_unit is None:
        length_unit = "code_length"
    if mass_unit is None:
        mass_unit = "code_mass"
    if time_unit is None:
        time_unit = "code_time"
    if velocity_unit is None:
        velocity_unit = "code_velocity"
    if magnetic_unit is None:
        magnetic_unit = "code_magnetic"

    particle_types = {}

    for grid in sfh.values():
        particle_types.update(set_particle_types(grid))

    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        parent_ids,
        number_of_particles,
        np.zeros(ngrids).reshape((ngrids, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "AMRGridData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = refine_by
    if np.all(domain_dimensions[1:] == 1):
        dimensionality = 1
    elif domain_dimensions[2] == 1:
        dimensionality = 2
    else:
        dimensionality = 3
    handler.dimensionality = dimensionality
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamDataset(
        handler,
        geometry=geometry,
        unit_system=unit_system,
        default_species_fields=default_species_fields,
    )
    return sds


def load_particles(
    data,
    length_unit=None,
    bbox=None,
    sim_time=None,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    geometry="cartesian",
    unit_system="cgs",
    data_source=None,
    default_species_fields=None,
):
    r"""Load a set of particles into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamParticleHandler`.

    This will allow a collection of particle data to be loaded directly into
    yt and analyzed as would any others.  This comes with several caveats:

    * There must be sufficient space in memory to contain all the particle
      data.
    * Parallelism will be disappointing or non-existent in most cases.
    * Fluid fields are not supported.

    Note: in order for the dataset to take advantage of SPH functionality,
    the following two fields must be provided:
    * ('io', 'density')
    * ('io', 'smoothing_length')

    Parameters
    ----------
    data : dict
        This is a dict of numpy arrays or (numpy array, unit name) tuples,
        where the keys are the field names. Particles positions must be named
        "particle_position_x", "particle_position_y", and "particle_position_z".
    length_unit : float
        Conversion factor from simulation length units to centimeters
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units of the length_unit
    sim_time : float, optional
        The simulation time in seconds
    mass_unit : float
        Conversion factor from simulation mass units to grams
    time_unit : float
        Conversion factor from simulation time units to seconds
    velocity_unit : float
        Conversion factor from simulation velocity units to cm/s
    magnetic_unit : float
        Conversion factor from simulation magnetic units to gauss
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    data_source : YTSelectionContainer, optional
        If set, parameters like `bbox`, `sim_time`, and code units are derived
        from it.
    default_species_fields : string, optional
        If set, default species fields are created for H and He which also
        determine the mean molecular weight. Options are "ionized" and "neutral".

    Examples
    --------

    >>> pos = [np.random.random(128 * 128 * 128) for i in range(3)]
    >>> data = dict(
    ...     particle_position_x=pos[0],
    ...     particle_position_y=pos[1],
    ...     particle_position_z=pos[2],
    ... )
    >>> bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    >>> ds = load_particles(data, 3.08e24, bbox=bbox)

    """
    from yt.frontends.stream.data_structures import (
        StreamDictFieldHandler,
        StreamHandler,
        StreamParticlesDataset,
    )
    from yt.frontends.stream.definitions import process_data, set_particle_types

    domain_dimensions = np.ones(3, "int32")
    nprocs = 1

    # Parse bounding box
    if data_source is not None:
        le, re = data_source.get_bbox()
        le = le.to_value("code_length")
        re = re.to_value("code_length")
        bbox = list(zip(le, re))
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], "float64")
    else:
        bbox = np.array(bbox)
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros(nprocs, dtype="int32").reshape((nprocs, 1))

    # Parse simulation time
    if data_source is not None:
        sim_time = data_source.ds.current_time
    if sim_time is None:
        sim_time = 0.0
    else:
        sim_time = float(sim_time)

    # Parse units
    def parse_unit(unit, dimension):
        if unit is None:
            unit = "code_" + dimension
            if data_source is not None:
                unit = getattr(data_source.ds, dimension + "_unit", unit)
        return unit

    length_unit = parse_unit(length_unit, "length")
    mass_unit = parse_unit(mass_unit, "mass")
    time_unit = parse_unit(time_unit, "time")
    velocity_unit = parse_unit(velocity_unit, "velocity")
    magnetic_unit = parse_unit(magnetic_unit, "magnetic")

    # Preprocess data
    field_units, data, _ = process_data(data)
    sfh = StreamDictFieldHandler()

    pdata = {}
    for key in data.keys():
        if not isinstance(key, tuple):
            field = ("io", key)
            mylog.debug("Reassigning '%s' to '%s'", key, field)
        else:
            field = key
        pdata[field] = data[key]
        sfh._additional_fields += (field,)
    data = pdata  # Drop reference count
    particle_types = set_particle_types(data)
    sfh.update({"stream_file": data})
    grid_left_edges = domain_left_edge
    grid_right_edges = domain_right_edge
    grid_dimensions = domain_dimensions.reshape(nprocs, 3).astype("int32")

    # I'm not sure we need any of this.
    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype="int64"),
        np.zeros(nprocs, dtype="int64").reshape(nprocs, 1),  # Temporary
        np.zeros(nprocs).reshape((nprocs, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "ParticleData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = 3
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamParticlesDataset(
        handler,
        geometry=geometry,
        unit_system=unit_system,
        default_species_fields=default_species_fields,
    )

    return sds


def load_hexahedral_mesh(
    data,
    connectivity,
    coordinates,
    length_unit=None,
    bbox=None,
    sim_time=0.0,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    geometry="cartesian",
    unit_system="cgs",
):
    r"""Load a hexahedral mesh of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.

    This should allow a semistructured grid of data to be loaded directly into
    yt and analyzed as would any others.  This comes with several caveats:

    * Units will be incorrect unless the data has already been converted to
      cgs.
    * Some functions may behave oddly, and parallelism will be
      disappointing or non-existent in most cases.
    * Particles may be difficult to integrate.

    Particle fields are detected as one-dimensional fields. The number of particles
    is set by the "number_of_particles" key in data.

    Parameters
    ----------
    data : dict
        This is a dict of numpy arrays, where the keys are the field names.
        There must only be one. Note that the data in the numpy arrays should
        define the cell-averaged value for of the quantity in in the hexahedral
        cell.
    connectivity : array_like
        This should be of size (N,8) where N is the number of zones.
    coordinates : array_like
        This should be of size (M,3) where M is the number of vertices
        indicated in the connectivity matrix.
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units of the length unit.
    sim_time : float, optional
        The simulation time in seconds
    mass_unit : string
        Unit to use for masses.  Defaults to unitless.
    time_unit : string
        Unit to use for times.  Defaults to unitless.
    velocity_unit : string
        Unit to use for velocities.  Defaults to unitless.
    magnetic_unit : string
        Unit to use for magnetic fields. Defaults to unitless.
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    geometry : string or tuple
        "cartesian", "cylindrical", "polar", "spherical", "geographic" or
        "spectral_cube".  Optionally, a tuple can be provided to specify the
        axis ordering -- for instance, to specify that the axis ordering should
        be z, x, y, this would be: ("cartesian", ("z", "x", "y")).  The same
        can be done for other coordinates, for instance:
        ("spherical", ("theta", "phi", "r")).

    """
    from yt.frontends.stream.data_structures import (
        StreamDictFieldHandler,
        StreamHandler,
        StreamHexahedralDataset,
    )
    from yt.frontends.stream.definitions import process_data, set_particle_types

    domain_dimensions = np.ones(3, "int32") * 2
    nprocs = 1
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], "float64")
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros(nprocs, dtype="int32").reshape((nprocs, 1))

    field_units, data, _ = process_data(data)
    sfh = StreamDictFieldHandler()

    particle_types = set_particle_types(data)

    sfh.update({"connectivity": connectivity, "coordinates": coordinates, 0: data})
    # Simple check for axis length correctness
    if len(data) > 0:
        fn = list(sorted(data))[0]
        array_values = data[fn]
        if array_values.size != connectivity.shape[0]:
            mylog.error(
                "Dimensions of array must be one fewer than the coordinate set."
            )
            raise RuntimeError
    grid_left_edges = domain_left_edge
    grid_right_edges = domain_right_edge
    grid_dimensions = domain_dimensions.reshape(nprocs, 3).astype("int32")

    if length_unit is None:
        length_unit = "code_length"
    if mass_unit is None:
        mass_unit = "code_mass"
    if time_unit is None:
        time_unit = "code_time"
    if velocity_unit is None:
        velocity_unit = "code_velocity"
    if magnetic_unit is None:
        magnetic_unit = "code_magnetic"

    # I'm not sure we need any of this.
    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype="int64"),
        np.zeros(nprocs, dtype="int64").reshape(nprocs, 1),  # Temporary
        np.zeros(nprocs).reshape((nprocs, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "HexahedralMeshData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = 3
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamHexahedralDataset(handler, geometry=geometry, unit_system=unit_system)

    return sds


def load_octree(
    octree_mask,
    data,
    bbox=None,
    sim_time=0.0,
    length_unit=None,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(True, True, True),
    over_refine_factor=1,
    partial_coverage=1,
    unit_system="cgs",
    default_species_fields=None,
):
    r"""Load an octree mask into yt.

    Octrees can be saved out by calling save_octree on an OctreeContainer.
    This enables them to be loaded back in.

    This will initialize an Octree of data.  Note that fluid fields will not
    work yet, or possibly ever.

    Parameters
    ----------
    octree_mask : np.ndarray[uint8_t]
        This is a depth-first refinement mask for an Octree.  It should be
        of size n_octs * 8 (but see note about the root oct below), where
        each item is 1 for an oct-cell being refined and 0 for it not being
        refined.  For over_refine_factors != 1, the children count will
        still be 8, so there will still be n_octs * 8 entries. Note that if
        the root oct is not refined, there will be only one entry
        for the root, so the size of the mask will be (n_octs - 1)*8 + 1.
    data : dict
        A dictionary of 1D arrays.  Note that these must of the size of the
        number of "False" values in the ``octree_mask``.
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units of length
    sim_time : float, optional
        The simulation time in seconds
    length_unit : string
        Unit to use for lengths.  Defaults to unitless.
    mass_unit : string
        Unit to use for masses.  Defaults to unitless.
    time_unit : string
        Unit to use for times.  Defaults to unitless.
    velocity_unit : string
        Unit to use for velocities.  Defaults to unitless.
    magnetic_unit : string
        Unit to use for magnetic fields. Defaults to unitless.
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    partial_coverage : boolean
        Whether or not an oct can be refined cell-by-cell, or whether all
        8 get refined.
    default_species_fields : string, optional
        If set, default species fields are created for H and He which also
        determine the mean molecular weight. Options are "ionized" and "neutral".

    Example
    -------

    >>> import numpy as np
    >>> oct_mask = np.zeros(25)
    ... oct_mask[[0,  5,  7, 16]] = 8
    >>> octree_mask = np.array(oct_mask, dtype=np.uint8)
    >>> quantities = {}
    >>> quantities["gas", "density"] = np.random.random((22, 1))
    >>> bbox = np.array([[-10.0, 10.0], [-10.0, 10.0], [-10.0, 10.0]])

    >>> ds = load_octree(
    ...     octree_mask=octree_mask,
    ...     data=quantities,
    ...     bbox=bbox,
    ...     over_refine_factor=0,
    ...     partial_coverage=0,
    ... )

    """
    from yt.frontends.stream.data_structures import (
        StreamDictFieldHandler,
        StreamHandler,
        StreamOctreeDataset,
    )
    from yt.frontends.stream.definitions import process_data, set_particle_types

    if not isinstance(octree_mask, np.ndarray) or octree_mask.dtype != np.uint8:
        raise TypeError("octree_mask should be a Numpy array with type uint8")

    nz = 1 << (over_refine_factor)
    domain_dimensions = np.array([nz, nz, nz])
    nprocs = 1
    if bbox is None:
        bbox = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], "float64")
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros(nprocs, dtype="int32").reshape((nprocs, 1))

    field_units, data, _ = process_data(data)
    sfh = StreamDictFieldHandler()

    particle_types = set_particle_types(data)

    sfh.update({0: data})
    grid_left_edges = domain_left_edge
    grid_right_edges = domain_right_edge
    grid_dimensions = domain_dimensions.reshape(nprocs, 3).astype("int32")

    if length_unit is None:
        length_unit = "code_length"
    if mass_unit is None:
        mass_unit = "code_mass"
    if time_unit is None:
        time_unit = "code_time"
    if velocity_unit is None:
        velocity_unit = "code_velocity"
    if magnetic_unit is None:
        magnetic_unit = "code_magnetic"

    # I'm not sure we need any of this.
    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype="int64"),
        np.zeros(nprocs, dtype="int64").reshape(nprocs, 1),  # Temporary
        np.zeros(nprocs).reshape((nprocs, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "OctreeData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = 3
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamOctreeDataset(
        handler, unit_system=unit_system, default_species_fields=default_species_fields
    )
    sds.octree_mask = octree_mask
    sds.partial_coverage = partial_coverage
    sds.over_refine_factor = over_refine_factor

    return sds


def load_unstructured_mesh(
    connectivity,
    coordinates,
    node_data=None,
    elem_data=None,
    length_unit=None,
    bbox=None,
    sim_time=0.0,
    mass_unit=None,
    time_unit=None,
    velocity_unit=None,
    magnetic_unit=None,
    periodicity=(False, False, False),
    geometry="cartesian",
    unit_system="cgs",
):
    r"""Load an unstructured mesh of data into yt as a
    :class:`~yt.frontends.stream.data_structures.StreamHandler`.

    This should allow an unstructured mesh data to be loaded directly into
    yt and analyzed as would any others.  Not all functionality for
    visualization will be present, and some analysis functions may not yet have
    been implemented.

    Particle fields are detected as one-dimensional fields. The number of
    particles is set by the "number_of_particles" key in data.

    In the parameter descriptions below, a "vertex" is a 3D point in space, an
    "element" is a single polyhedron whose location is defined by a set of
    vertices, and a "mesh" is a set of polyhedral elements, each with the same
    number of vertices.

    Parameters
    ----------

    connectivity : list of array_like or array_like
        This should either be a single 2D array or list of 2D arrays.  If this
        is a list, each element in the list corresponds to the connectivity
        information for a distinct mesh. Each array can have different
        connectivity length and should be of shape (N,M) where N is the number
        of elements and M is the number of vertices per element.
    coordinates : array_like
        The 3D coordinates of mesh vertices. This should be of size (L, D) where
        L is the number of vertices and D is the number of coordinates per vertex
        (the spatial dimensions of the dataset). Currently this must be either 2 or 3.
        When loading more than one mesh, the data for each mesh should be concatenated
        into a single coordinates array.
    node_data : dict or list of dicts
        For a single mesh, a dict mapping field names to 2D numpy arrays,
        representing data defined at element vertices. For multiple meshes,
        this must be a list of dicts.  Note that these are not the values as a
        function of the coordinates, but of the connectivity.  Their shape
        should be the same as the connectivity.  This means that if the data is
        in the shape of the coordinates, you may need to reshape them using the
        `connectivity` array as an index.
    elem_data : dict or list of dicts
        For a single mesh, a dict mapping field names to 1D numpy arrays, where
        each array has a length equal to the number of elements. The data
        must be defined at the center of each mesh element and there must be
        only one data value for each element. For multiple meshes, this must be
        a list of dicts, with one dict for each mesh.
    bbox : array_like (xdim:zdim, LE:RE), optional
        Size of computational domain in units of the length unit.
    sim_time : float, optional
        The simulation time in seconds
    length_unit : string
        Unit to use for length.  Defaults to unitless.
    mass_unit : string
        Unit to use for masses.  Defaults to unitless.
    time_unit : string
        Unit to use for times.  Defaults to unitless.
    velocity_unit : string
        Unit to use for velocities.  Defaults to unitless.
    magnetic_unit : string
        Unit to use for magnetic fields. Defaults to unitless.
    periodicity : tuple of booleans
        Determines whether the data will be treated as periodic along
        each axis
    geometry : string or tuple
        "cartesian", "cylindrical", "polar", "spherical", "geographic" or
        "spectral_cube".  Optionally, a tuple can be provided to specify the
        axis ordering -- for instance, to specify that the axis ordering should
        be z, x, y, this would be: ("cartesian", ("z", "x", "y")).  The same
        can be done for other coordinates, for instance:
        ("spherical", ("theta", "phi", "r")).

    Examples
    --------

    Load a simple mesh consisting of two tets.

      >>> # Coordinates for vertices of two tetrahedra
      >>> coordinates = np.array(
      ...     [
      ...         [0.0, 0.0, 0.5],
      ...         [0.0, 1.0, 0.5],
      ...         [0.5, 1, 0.5],
      ...         [0.5, 0.5, 0.0],
      ...         [0.5, 0.5, 1.0],
      ...     ]
      ... )
      >>> # The indices in the coordinates array of mesh vertices.
      >>> # This mesh has two elements.
      >>> connectivity = np.array([[0, 1, 2, 4], [0, 1, 2, 3]])

      >>> # Field data defined at the centers of the two mesh elements.
      >>> elem_data = {("connect1", "elem_field"): np.array([1, 2])}

      >>> # Field data defined at node vertices
      >>> node_data = {
      ...     ("connect1", "node_field"): np.array(
      ...         [[0.0, 1.0, 2.0, 4.0], [0.0, 1.0, 2.0, 3.0]]
      ...     )
      ... }

      >>> ds = load_unstructured_mesh(
      ...     connectivity, coordinates, elem_data=elem_data, node_data=node_data
      ... )
    """
    from yt.frontends.exodus_ii.util import get_num_pseudo_dims
    from yt.frontends.stream.data_structures import (
        StreamDictFieldHandler,
        StreamHandler,
        StreamUnstructuredMeshDataset,
    )
    from yt.frontends.stream.definitions import process_data, set_particle_types

    dimensionality = coordinates.shape[1]
    domain_dimensions = np.ones(3, "int32") * 2
    nprocs = 1

    if elem_data is None and node_data is None:
        raise RuntimeError("No data supplied in load_unstructured_mesh.")

    connectivity = list(always_iterable(connectivity, base_type=np.ndarray))
    num_meshes = max(1, len(connectivity))

    elem_data = list(always_iterable(elem_data, base_type=dict)) or [{}] * num_meshes
    node_data = list(always_iterable(node_data, base_type=dict)) or [{}] * num_meshes

    data = [{} for i in range(num_meshes)]
    for elem_dict, data_dict in zip(elem_data, data):
        for field, values in elem_dict.items():
            data_dict[field] = values
    for node_dict, data_dict in zip(node_data, data):
        for field, values in node_dict.items():
            data_dict[field] = values

    if bbox is None:
        bbox = [
            [
                coordinates[:, i].min() - 0.1 * abs(coordinates[:, i].min()),
                coordinates[:, i].max() + 0.1 * abs(coordinates[:, i].max()),
            ]
            for i in range(dimensionality)
        ]

    if dimensionality < 3:
        bbox.append([0.0, 1.0])
    if dimensionality < 2:
        bbox.append([0.0, 1.0])

    # handle pseudo-dims here
    num_pseudo_dims = get_num_pseudo_dims(coordinates)
    dimensionality -= num_pseudo_dims
    for i in range(dimensionality, 3):
        bbox[i][0] = 0.0
        bbox[i][1] = 1.0

    bbox = np.array(bbox, dtype=np.float64)
    domain_left_edge = np.array(bbox[:, 0], "float64")
    domain_right_edge = np.array(bbox[:, 1], "float64")
    grid_levels = np.zeros(nprocs, dtype="int32").reshape((nprocs, 1))

    field_units = {}
    particle_types = {}
    sfh = StreamDictFieldHandler()

    sfh.update({"connectivity": connectivity, "coordinates": coordinates})
    for i, d in enumerate(data):
        _f_unit, _data, _ = process_data(d)
        field_units.update(_f_unit)
        sfh[i] = _data
        particle_types.update(set_particle_types(d))

    grid_left_edges = domain_left_edge
    grid_right_edges = domain_right_edge
    grid_dimensions = domain_dimensions.reshape(nprocs, 3).astype("int32")

    if length_unit is None:
        length_unit = "code_length"
    if mass_unit is None:
        mass_unit = "code_mass"
    if time_unit is None:
        time_unit = "code_time"
    if velocity_unit is None:
        velocity_unit = "code_velocity"
    if magnetic_unit is None:
        magnetic_unit = "code_magnetic"

    # I'm not sure we need any of this.
    handler = StreamHandler(
        grid_left_edges,
        grid_right_edges,
        grid_dimensions,
        grid_levels,
        -np.ones(nprocs, dtype="int64"),
        np.zeros(nprocs, dtype="int64").reshape(nprocs, 1),  # Temporary
        np.zeros(nprocs).reshape((nprocs, 1)),
        sfh,
        field_units,
        (length_unit, mass_unit, time_unit, velocity_unit, magnetic_unit),
        particle_types=particle_types,
        periodicity=periodicity,
    )

    handler.name = "UnstructuredMeshData"
    handler.domain_left_edge = domain_left_edge
    handler.domain_right_edge = domain_right_edge
    handler.refine_by = 2
    handler.dimensionality = dimensionality
    handler.domain_dimensions = domain_dimensions
    handler.simulation_time = sim_time
    handler.cosmology_simulation = 0

    sds = StreamUnstructuredMeshDataset(
        handler, geometry=geometry, unit_system=unit_system
    )

    fluid_types = ["all"]
    for i in range(1, num_meshes + 1):
        fluid_types += ["connect%d" % i]
    sds.fluid_types = tuple(fluid_types)

    def flatten(l):
        return [item for sublist in l for item in sublist]

    sds._node_fields = flatten([[f[1] for f in m] for m in node_data if m])
    sds._elem_fields = flatten([[f[1] for f in m] for m in elem_data if m])
    sds.default_field = [f for f in sds.field_list if f[0] == "connect1"][-1]
    sds.default_fluid_type = sds.default_field[0]
    return sds


# --- Loader for yt sample datasets ---
def load_sample(
    fn: Optional[str] = None, *, progressbar: bool = True, timeout=None, **kwargs
):
    r"""
    Load sample data with yt.

    This is a simple wrapper around :func:`~yt.loaders.load` to include fetching
    data with pooch from remote source.

    The data registry table can be retrieved and visualized using
    :func:`~yt.sample_data.api.get_data_registry_table`.
    The `filename` column contains usable keys that can be passed
    as the first positional argument to load_sample.
    Some data samples contain series of datasets. It may be required to
    supply the relative path to a specific dataset.

    Parameters
    ----------

    fn: str
        The `filename` of the dataset to load, as defined in the data registry
        table.

    progressbar: bool
        display a progress bar (tqdm).

    timeout: float or int (optional)
        Maximal waiting time, in seconds, after which download is aborted.
        `None` means "no limit". This parameter is directly passed to down to
        requests.get via pooch.HTTPDownloader

    Notes
    -----

    - This function is experimental as of yt 4.0.0, do not rely on its exact behaviour.
    - Any additional keyword argument is passed down to :func:`~yt.loaders.load`.
    - In case of collision with predefined keyword arguments as set in
      the data registry, the ones passed to this function take priority.
    - Datasets with slashes '/' in their names can safely be used even on Windows.
      On the contrary, paths using backslashes '\' won't work outside of Windows, so
      it is recommended to favour the UNIX convention ('/') in scripts that are meant
      to be cross-platform.
    - This function requires pandas and pooch.
    - Corresponding sample data live at https://yt-project.org/data

    """

    if fn is None:
        print(
            "One can see which sample datasets are available at: https://yt-project.org/data\n"
            "or alternatively by running: yt.sample_data.api.get_data_registry_table()",
            file=sys.stderr,
        )
        return None

    from yt.sample_data.api import (
        _download_sample_data_file,
        _get_test_data_dir_path,
        get_data_registry_table,
    )

    pooch_logger = pooch.utils.get_logger()

    # normalize path for platform portability
    # for consistency with yt.load, we also convert to str explicitly,
    # which gives us support Path objects for free
    fn = str(fn).replace("/", os.path.sep)

    topdir, _, specific_file = fn.partition(os.path.sep)

    registry_table = get_data_registry_table()

    known_names: List[str] = registry_table.dropna()["filename"].to_list()
    if topdir not in known_names:
        msg = f"'{topdir}' is not an available dataset."
        lexical_distances: List[Tuple[str, int]] = [
            (name, levenshtein_distance(name, topdir)) for name in known_names
        ]
        suggestions: List[str] = [name for name, dist in lexical_distances if dist < 4]
        if len(suggestions) == 1:
            msg += f" Did you mean '{suggestions[0]}' ?"
        elif suggestions:
            msg += " Did you mean to type any of the following ?\n\n    "
            msg += "\n    ".join(f"'{_}'" for _ in suggestions)
        raise ValueError(msg)

    # PR 3089
    # note: in the future the registry table should be reindexed
    # so that the following line can be replaced with
    #
    # specs = registry_table.loc[fn]
    #
    # however we don't want to do it right now because the "filename" column is
    # currently incomplete
    specs = registry_table.query(f"`filename` == '{topdir}'").iloc[0]

    load_name = specific_file or specs["load_name"] or ""

    if not isinstance(specs["load_kwargs"], dict):
        raise ValueError(
            "The requested dataset seems to be improperly registered.\n"
            "Tip: the entry in yt/sample_data_registry.json may be inconsistent with "
            "https://github.com/yt-project/website/blob/master/data/datafiles.json\n"
            "Please report this to https://github.com/yt-project/yt/issues/new"
        )

    kwargs = {**specs["load_kwargs"], **kwargs}

    save_dir = _get_test_data_dir_path()

    data_path = save_dir.joinpath(fn)
    if save_dir.joinpath(topdir).exists():
        # if the data is already available locally, `load_sample`
        # only acts as a thin wrapper around `load`
        if load_name and os.sep not in fn:
            data_path = data_path.joinpath(load_name)
        mylog.info("Sample dataset found in '%s'", data_path)
        if timeout is not None:
            mylog.info("Ignoring the `timeout` keyword argument received.")
        return load(data_path, **kwargs)

    mylog.info("'%s' is not available locally. Looking up online.", fn)

    # effectively silence the pooch's logger and create our own log instead
    pooch_logger.setLevel(100)
    mylog.info("Downloading from %s", specs["url"])

    # downloading via a pooch.Pooch instance behind the scenes
    filename = urlsplit(specs["url"]).path.split("/")[-1]

    tmp_file = _download_sample_data_file(
        filename, progressbar=progressbar, timeout=timeout
    )

    # pooch has functionalities to unpack downloaded archive files,
    # but it needs to be told in advance that we are downloading a tarball.
    # Since that information is not necessarily trivial to guess from the filename,
    # we rely on the standard library to perform a conditional unpacking instead.
    if tarfile.is_tarfile(tmp_file):
        mylog.info("Untaring downloaded file to '%s'", save_dir)
        with tarfile.open(tmp_file) as fh:
            fh.extractall(save_dir)
        os.remove(tmp_file)
    else:
        os.replace(tmp_file, save_dir)

    loadable_path = Path.joinpath(save_dir, fn)
    if load_name not in str(loadable_path):
        loadable_path = loadable_path.joinpath(load_name, specific_file)

    return load(loadable_path, **kwargs)

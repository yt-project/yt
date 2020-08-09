from collections import defaultdict

import numpy as np

from yt.funcs import iterable
from yt.geometry.grid_container import GridTree, MatchPointsToGrids
from yt.units.yt_array import uconcatenate
from yt.utilities.exceptions import (
    YTInconsistentGridFieldShape,
    YTInconsistentGridFieldShapeGridDims,
    YTInconsistentParticleFieldShape,
)
from yt.utilities.flagging_methods import FlaggingGrid
from yt.utilities.logger import ytLogger as mylog

from .fields import StreamFieldInfo


def assign_particle_data(ds, pdata, bbox):

    """
    Assign particle data to the grids using MatchPointsToGrids. This
    will overwrite any existing particle data, so be careful!
    """

    for ptype in ds.particle_types_raw:
        check_fields = [(ptype, "particle_position_x"), (ptype, "particle_position")]
        if all(f not in pdata for f in check_fields):
            pdata_ftype = {}
            for f in [k for k in sorted(pdata)]:
                if not hasattr(pdata[f], "shape"):
                    continue
                if f == "number_of_particles":
                    continue
                mylog.debug("Reassigning '%s' to ('%s','%s')", f, ptype, f)
                pdata_ftype[ptype, f] = pdata.pop(f)
            pdata_ftype.update(pdata)
            pdata = pdata_ftype

    # Note: what we need to do here is a bit tricky.  Because occasionally this
    # gets called before we property handle the field detection, we cannot use
    # any information about the index.  Fortunately for us, we can generate
    # most of the GridTree utilizing information we already have from the
    # stream handler.

    if len(ds.stream_handler.fields) > 1:
        pdata.pop("number_of_particles", None)
        num_grids = len(ds.stream_handler.fields)
        parent_ids = ds.stream_handler.parent_ids
        num_children = np.zeros(num_grids, dtype="int64")
        # We're going to do this the slow way
        mask = np.empty(num_grids, dtype="bool")
        for i in range(num_grids):
            np.equal(parent_ids, i, mask)
            num_children[i] = mask.sum()
        levels = ds.stream_handler.levels.astype("int64").ravel()
        grid_tree = GridTree(
            num_grids,
            ds.stream_handler.left_edges,
            ds.stream_handler.right_edges,
            ds.stream_handler.dimensions,
            ds.stream_handler.parent_ids,
            levels,
            num_children,
        )

        grid_pdata = []
        for _ in range(num_grids):
            grid = {"number_of_particles": 0}
            grid_pdata.append(grid)

        for ptype in ds.particle_types_raw:
            if (ptype, "particle_position_x") in pdata:
                x, y, z = (pdata[ptype, "particle_position_%s" % ax] for ax in "xyz")
            elif (ptype, "particle_position") in pdata:
                x, y, z = pdata[ptype, "particle_position"].T
            else:
                raise KeyError(
                    "Cannot decompose particle data without position fields!"
                )
            pts = MatchPointsToGrids(grid_tree, len(x), x, y, z)
            particle_grid_inds = pts.find_points_in_tree()
            (assigned_particles,) = (particle_grid_inds >= 0).nonzero()
            num_particles = particle_grid_inds.size
            num_unassigned = num_particles - assigned_particles.size
            if num_unassigned > 0:
                eps = np.finfo(x.dtype).eps
                s = np.array(
                    [
                        [x.min() - eps, x.max() + eps],
                        [y.min() - eps, y.max() + eps],
                        [z.min() - eps, z.max() + eps],
                    ]
                )
                sug_bbox = [
                    [min(bbox[0, 0], s[0, 0]), max(bbox[0, 1], s[0, 1])],
                    [min(bbox[1, 0], s[1, 0]), max(bbox[1, 1], s[1, 1])],
                    [min(bbox[2, 0], s[2, 0]), max(bbox[2, 1], s[2, 1])],
                ]
                mylog.warning(
                    "Discarding %s particles (out of %s) that are outside "
                    "bounding box. Set bbox=%s to avoid this in the future.",
                    num_unassigned,
                    num_particles,
                    sug_bbox,
                )
                particle_grid_inds = particle_grid_inds[assigned_particles]
                x = x[assigned_particles]
                y = y[assigned_particles]
                z = z[assigned_particles]
            idxs = np.argsort(particle_grid_inds)
            particle_grid_count = np.bincount(
                particle_grid_inds.astype("intp"), minlength=num_grids
            )
            particle_indices = np.zeros(num_grids + 1, dtype="int64")
            if num_grids > 1:
                np.add.accumulate(
                    particle_grid_count.squeeze(), out=particle_indices[1:]
                )
            else:
                particle_indices[1] = particle_grid_count.squeeze()
            for i, pcount in enumerate(particle_grid_count):
                grid_pdata[i]["number_of_particles"] += pcount
                start = particle_indices[i]
                end = particle_indices[i + 1]
                for key in pdata.keys():
                    if key[0] == ptype:
                        grid_pdata[i][key] = pdata[key][idxs][start:end]

    else:
        grid_pdata = [pdata]

    for pd, gi in zip(grid_pdata, sorted(ds.stream_handler.fields)):
        ds.stream_handler.fields[gi].update(pd)
        ds.stream_handler.particle_types.update(set_particle_types(pd))
        npart = ds.stream_handler.fields[gi].pop("number_of_particles", 0)
        ds.stream_handler.particle_count[gi] = npart


def process_data(data, grid_dims=None):
    new_data, field_units = {}, {}
    for field, val in data.items():
        # val is a data array
        if isinstance(val, np.ndarray):
            # val is a YTArray
            if hasattr(val, "units"):
                field_units[field] = val.units
                new_data[field] = val.copy().d
            # val is a numpy array
            else:
                field_units[field] = ""
                new_data[field] = val.copy()

        # val is a tuple of (data, units)
        elif isinstance(val, tuple) and len(val) == 2:
            try:
                assert isinstance(field, (str, tuple)), "Field name is not a string!"
                assert isinstance(val[0], np.ndarray), "Field data is not an ndarray!"
                assert isinstance(val[1], str), "Unit specification is not a string!"
                field_units[field] = val[1]
                new_data[field] = val[0]
            except AssertionError as e:
                raise RuntimeError("The data dict appears to be invalid.\n" + str(e))

        # val is a list of data to be turned into an array
        elif iterable(val):
            field_units[field] = ""
            new_data[field] = np.asarray(val)

        else:
            raise RuntimeError(
                "The data dict appears to be invalid. "
                "The data dictionary must map from field "
                "names to (numpy array, unit spec) tuples. "
            )

    data = new_data

    # At this point, we have arrays for all our fields
    new_data = {}
    for field in data:
        n_shape = len(data[field].shape)
        if isinstance(field, tuple):
            new_field = field
        elif n_shape in (1, 2):
            new_field = ("io", field)
        elif n_shape == 3:
            new_field = ("stream", field)
        else:
            raise RuntimeError
        new_data[new_field] = data[field]
        field_units[new_field] = field_units.pop(field)
        known_fields = (
            StreamFieldInfo.known_particle_fields + StreamFieldInfo.known_other_fields
        )
        # We do not want to override any of the known ones, if it's not
        # overridden here.
        if (
            any(f[0] == new_field[1] for f in known_fields)
            and field_units[new_field] == ""
        ):
            field_units.pop(new_field)
    data = new_data
    # Sanity checking that all fields have the same dimensions.
    g_shapes = []
    p_shapes = defaultdict(list)
    for field in data:
        f_shape = data[field].shape
        n_shape = len(f_shape)
        if n_shape in (1, 2):
            p_shapes[field[0]].append((field[1], f_shape[0]))
        elif n_shape == 3:
            g_shapes.append((field, f_shape))
    if len(g_shapes) > 0:
        g_s = np.array([s[1] for s in g_shapes])
        if not np.all(g_s == g_s[0]):
            raise YTInconsistentGridFieldShape(g_shapes)
        if grid_dims is not None:
            if not np.all(g_s == grid_dims):
                raise YTInconsistentGridFieldShapeGridDims(g_shapes, grid_dims)
    if len(p_shapes) > 0:
        for ptype, p_shape in p_shapes.items():
            p_s = np.array([s[1] for s in p_shape])
            if not np.all(p_s == p_s[0]):
                raise YTInconsistentParticleFieldShape(ptype, p_shape)
    # Now that we know the particle fields are consistent, determine the number
    # of particles.
    if len(p_shapes) > 0:
        number_of_particles = np.sum([s[0][1] for s in p_shapes.values()])
    else:
        number_of_particles = 0
    return field_units, data, number_of_particles


def refine_amr(base_ds, refinement_criteria, fluid_operators, max_level, callback=None):
    r"""Given a base dataset, repeatedly apply refinement criteria and
    fluid operators until a maximum level is reached.

    Parameters
    ----------
    base_ds : ~yt.data_objects.static_output.Dataset
        This is any static output.  It can also be a stream static output, for
        instance as returned by `yt.loaders.load_uniform_data`.
    refinement_critera : list of :class:`~yt.utilities.flagging_methods.FlaggingMethod`
        These criteria will be applied in sequence to identify cells that need
        to be refined.
    fluid_operators : list of :class:`~yt.utilities.initial_conditions.FluidOperator`
        These fluid operators will be applied in sequence to all resulting
        grids.
    max_level : int
        The maximum level to which the data will be refined
    callback : function, optional
        A function that will be called at the beginning of each refinement
        cycle, with the current dataset.

    Examples
    --------
    >>> domain_dims = (32, 32, 32)
    >>> data = np.zeros(domain_dims) + 0.25
    >>> fo = [ic.CoredSphere(0.05, 0.3, [0.7,0.4,0.75], {"Density": (0.25, 100.0)})]
    >>> rc = [fm.flagging_method_registry["overdensity"](8.0)]
    >>> ug = load_uniform_grid({'Density': data}, domain_dims, 1.0)
    >>> ds = refine_amr(ug, rc, fo, 5)
    """
    from .loaders import load_amr_grids

    # If we have particle data, set it aside for now

    number_of_particles = np.sum(
        [grid.NumberOfParticles for grid in base_ds.index.grids]
    )

    if number_of_particles > 0:
        pdata = {}
        for field in base_ds.field_list:
            if not isinstance(field, tuple):
                field = ("unknown", field)
            fi = base_ds._get_field_info(*field)
            if (
                fi.sampling_type == "particle"
                and field[0] in base_ds.particle_types_raw
            ):
                pdata[field] = uconcatenate(
                    [grid[field] for grid in base_ds.index.grids]
                )
        pdata["number_of_particles"] = number_of_particles

    last_gc = base_ds.index.num_grids
    cur_gc = -1
    ds = base_ds
    bbox = np.array(
        [(ds.domain_left_edge[i], ds.domain_right_edge[i]) for i in range(3)]
    )
    while ds.index.max_level < max_level and last_gc != cur_gc:
        mylog.info("Refining another level.  Current max level: %s", ds.index.max_level)
        last_gc = ds.index.grids.size
        for m in fluid_operators:
            m.apply(ds)
        if callback is not None:
            callback(ds)
        grid_data = []
        for g in ds.index.grids:
            gd = dict(
                left_edge=g.LeftEdge,
                right_edge=g.RightEdge,
                level=g.Level,
                dimensions=g.ActiveDimensions,
            )
            for field in ds.field_list:
                if not isinstance(field, tuple):
                    field = ("unknown", field)
                fi = ds._get_field_info(*field)
                if not fi.sampling_type == "particle":
                    gd[field] = g[field]
            grid_data.append(gd)
            if g.Level < ds.index.max_level:
                continue
            fg = FlaggingGrid(g, refinement_criteria)
            nsg = fg.find_subgrids()
            for sg in nsg:
                LE = sg.left_index * g.dds + ds.domain_left_edge
                dims = sg.dimensions * ds.refine_by
                grid = ds.smoothed_covering_grid(g.Level + 1, LE, dims)
                gd = dict(
                    left_edge=LE,
                    right_edge=grid.right_edge,
                    level=g.Level + 1,
                    dimensions=dims,
                )
                for field in ds.field_list:
                    if not isinstance(field, tuple):
                        field = ("unknown", field)
                    fi = ds._get_field_info(*field)
                    if not fi.sampling_type == "particle":
                        gd[field] = grid[field]
                grid_data.append(gd)

        ds = load_amr_grids(grid_data, ds.domain_dimensions, bbox=bbox)

        ds.particle_types_raw = base_ds.particle_types_raw
        ds.particle_types = ds.particle_types_raw

        # Now figure out where the particles go
        if number_of_particles > 0:
            # This will update the stream handler too
            assign_particle_data(ds, pdata, bbox)

        cur_gc = ds.index.num_grids

    return ds


def set_particle_types(data):
    particle_types = {}
    for key in data.keys():
        if key == "number_of_particles":
            continue
        if len(data[key].shape) == 1:
            particle_types[key] = True
        else:
            particle_types[key] = False
    return particle_types

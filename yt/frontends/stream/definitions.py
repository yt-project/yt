from collections import defaultdict

import numpy as np

from yt.funcs import is_sequence
from yt.geometry.grid_container import GridTree, MatchPointsToGrids
from yt.utilities.exceptions import (
    YTInconsistentGridFieldShape,
    YTInconsistentGridFieldShapeGridDims,
    YTInconsistentParticleFieldShape,
)
from yt.utilities.logger import ytLogger as mylog

from .fields import StreamFieldInfo


def assign_particle_data(ds, pdata, bbox):

    """
    Assign particle data to the grids using MatchPointsToGrids. This
    will overwrite any existing particle data, so be careful!
    """

    particle_index_fields = [
        f"particle_position_{ax}" for ax in ds.coordinates.axis_order
    ]
    for ptype in ds.particle_types_raw:
        check_fields = [(ptype, pi_field) for pi_field in particle_index_fields]
        check_fields.append((ptype, "particle_position"))
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
            particle_index_fields = [
                f"particle_position_{ax}" for ax in ds.coordinates.axis_order
            ]

        for ptype in ds.particle_types_raw:
            if (ptype, "particle_position_x") in pdata:
                # we call them x, y, z even though they may be different field names
                x, y, z = (pdata[ptype, pi_field] for pi_field in particle_index_fields)
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


def process_data(data, grid_dims=None, allow_callables=True):
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
                valid_data = isinstance(val[0], np.ndarray)
                if allow_callables:
                    valid_data = valid_data or callable(val[0])
                assert isinstance(field, (str, tuple)), "Field name is not a string!"
                assert (
                    valid_data
                ), "Field data is not an ndarray or callable (with nproc == 1)!"
                assert isinstance(val[1], str), "Unit specification is not a string!"
                field_units[field] = val[1]
                new_data[field] = val[0]
            except AssertionError as e:
                raise RuntimeError("The data dict appears to be invalid.\n" + str(e))

        # val is a list of data to be turned into an array
        elif is_sequence(val):
            field_units[field] = ""
            new_data[field] = np.asarray(val)

        elif callable(val):
            if not allow_callables:
                raise RuntimeError(
                    "Callable functions can not be specified "
                    "in conjunction with nprocs > 1."
                )
            field_units[field] = ""
            new_data[field] = val
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
        n_shape = 3
        if not callable(data[field]):
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
        if callable(data[field]):
            continue
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


def set_particle_types(data):
    particle_types = {}
    for key in data.keys():
        if key == "number_of_particles":
            continue
        elif callable(data[key]):
            particle_types[key] = False
        elif len(data[key].shape) == 1:
            particle_types[key] = True
        else:
            particle_types[key] = False
    return particle_types

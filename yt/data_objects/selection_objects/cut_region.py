import re

import numpy as np
from more_itertools import always_iterable

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.funcs import iter_fields, validate_object, validate_sequence
from yt.geometry.selection_routines import points_in_cells
from yt.utilities.exceptions import YTIllDefinedCutRegion
from yt.utilities.on_demand_imports import _scipy


class YTCutRegion(YTSelectionContainer3D):
    """
    This is a data object designed to allow individuals to apply logical
    operations to fields and filter as a result of those cuts.

    Parameters
    ----------
    data_source : YTSelectionContainer3D
        The object to which cuts will be applied.
    conditionals : list of strings
        A list of conditionals that will be evaluated.  In the namespace
        available, these conditionals will have access to 'obj' which is a data
        object of unknown shape, and they must generate a boolean array.  For
        instance, conditionals = ["obj[('gas', 'temperature')] < 1e3"]

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("RedshiftOutput0005")
    >>> sp = ds.sphere("max", (1.0, "Mpc"))
    >>> cr = ds.cut_region(sp, ["obj[('gas', 'temperature')] < 1e3"])
    """

    _type_name = "cut_region"
    _con_args = ("base_object", "conditionals")
    _derived_quantity_chunking = "all"

    def __init__(
        self,
        data_source,
        conditionals,
        ds=None,
        field_parameters=None,
        base_object=None,
        locals=None,
    ):
        if locals is None:
            locals = {}
        validate_object(data_source, YTSelectionContainer)
        validate_sequence(conditionals)
        for condition in conditionals:
            validate_object(condition, str)
        validate_object(ds, Dataset)
        validate_object(field_parameters, dict)
        validate_object(base_object, YTSelectionContainer)

        self.conditionals = list(always_iterable(conditionals))
        if isinstance(data_source, YTCutRegion):
            # If the source is also a cut region, add its conditionals
            # and set the source to be its source.
            # Preserve order of conditionals.
            self.conditionals = data_source.conditionals + self.conditionals
            data_source = data_source.base_object

        super().__init__(
            data_source.center, ds, field_parameters, data_source=data_source
        )
        self.filter_fields = self._check_filter_fields()
        self.base_object = data_source
        self.locals = locals
        self._selector = None
        # Need to interpose for __getitem__, fwidth, fcoords, icoords, iwidth,
        # ires and get_data

    def _check_filter_fields(self):
        fields = []
        for cond in self.conditionals:
            for field in re.findall(r"\[([A-Za-z0-9_,.'\"\(\)]+)\]", cond):
                fd = field.replace('"', "").replace("'", "")
                if "," in fd:
                    fd = tuple(fd.strip("()").split(","))
                fd = self.ds._get_field_info(fd)
                if fd.sampling_type == "particle" or fd.is_sph_field:
                    raise RuntimeError(
                        f"cut_region requires a mesh-based field, "
                        f"but {fd.name} is a particle field! Use "
                        f"a particle filter instead. "
                    )
                fields.append(fd.name)
        return fields

    def chunks(self, fields, chunking_style, **kwargs):
        # We actually want to chunk the sub-chunk, not ourselves.  We have no
        # chunks to speak of, as we do not data IO.
        for chunk in self.index._chunk(self.base_object, chunking_style, **kwargs):
            with self.base_object._chunked_read(chunk):
                with self._chunked_read(chunk):
                    self.get_data(fields)
                    yield self

    def get_data(self, fields=None):
        fields = list(iter_fields(fields))
        self.base_object.get_data(fields)
        ind = self._cond_ind
        for field in fields:
            f = self.base_object[field]
            if f.shape != ind.shape:
                parent = getattr(self, "parent", self.base_object)
                self.field_data[field] = parent[field][self._part_ind(field[0])]
            else:
                self.field_data[field] = self.base_object[field][ind]

    @property
    def blocks(self):
        # We have to take a slightly different approach here.  Note that all
        # that .blocks has to yield is a 3D array and a mask.
        for obj, m in self.base_object.blocks:
            m = m.copy()
            with obj._field_parameter_state(self.field_parameters):
                for cond in self.conditionals:
                    ss = eval(cond)
                    m = np.logical_and(m, ss, m)
            if not np.any(m):
                continue
            yield obj, m

    @property
    def _cond_ind(self):
        ind = None
        obj = self.base_object
        locals = self.locals.copy()
        if "obj" in locals:
            raise RuntimeError(
                '"obj" has been defined in the "locals" ; '
                "this is not supported, please rename the variable."
            )
        locals["obj"] = obj
        with obj._field_parameter_state(self.field_parameters):
            for cond in self.conditionals:
                res = eval(cond, locals)
                if ind is None:
                    ind = res
                if ind.shape != res.shape:
                    raise YTIllDefinedCutRegion(self.conditionals)
                np.logical_and(res, ind, ind)
        return ind

    def _part_ind_KDTree(self, ptype):
        """Find the particles in cells using a KDTree approach."""
        parent = getattr(self, "parent", self.base_object)
        units = "code_length"

        pos = np.stack(
            [
                self[("index", "x")].to(units),
                self[("index", "y")].to(units),
                self[("index", "z")].to(units),
            ],
            axis=1,
        ).value
        dx = np.stack(
            [
                self[("index", "dx")].to(units),
                self[("index", "dy")].to(units),
                self[("index", "dz")].to(units),
            ],
            axis=1,
        ).value
        ppos = np.stack(
            [
                parent[(ptype, "particle_position_x")],
                parent[(ptype, "particle_position_y")],
                parent[(ptype, "particle_position_z")],
            ],
            axis=1,
        ).value

        mask = np.zeros(ppos.shape[0], dtype=bool)
        levels = self[("index", "grid_level")].astype("int32").value
        if levels.size == 0:
            return mask

        levelmin = levels.min()
        levelmax = levels.max()

        for lvl in range(levelmax, levelmin - 1, -1):
            # Filter out cells not in the current level
            lvl_mask = levels == lvl
            dx_loc = dx[lvl_mask]
            pos_loc = pos[lvl_mask]

            grid_tree = _scipy.spatial.cKDTree(pos_loc, boxsize=1)

            # Compute closest cell for all remaining particles
            dist, icell = grid_tree.query(
                ppos[~mask], distance_upper_bound=dx_loc.max(), p=np.inf
            )
            mask_loc = np.isfinite(dist[:])

            # Check that particles within dx of a cell are in it
            i = icell[mask_loc]
            dist = np.abs(ppos[~mask][mask_loc, :] - pos_loc[i])
            tmp_mask = np.all(dist <= (dx_loc[i] / 2), axis=1)

            mask_loc[mask_loc] = tmp_mask

            # Update the particle mask with particles found at this level
            mask[~mask] |= mask_loc

        return mask

    def _part_ind_brute_force(self, ptype):
        parent = getattr(self, "parent", self.base_object)
        units = "code_length"
        mask = points_in_cells(
            self[("index", "x")].to(units),
            self[("index", "y")].to(units),
            self[("index", "z")].to(units),
            self[("index", "dx")].to(units),
            self[("index", "dy")].to(units),
            self[("index", "dz")].to(units),
            parent[(ptype, "particle_position_x")].to(units),
            parent[(ptype, "particle_position_y")].to(units),
            parent[(ptype, "particle_position_z")].to(units),
        )

        return mask

    def _part_ind(self, ptype):
        # If scipy is installed, use the fast KD tree
        # implementation. Else, fall back onto the direct
        # brute-force algorithm.
        try:
            _scipy.spatial.KDTree
            return self._part_ind_KDTree(ptype)
        except ImportError:
            return self._part_ind_brute_force(ptype)

    @property
    def icoords(self):
        return self.base_object.icoords[self._cond_ind, :]

    @property
    def fcoords(self):
        return self.base_object.fcoords[self._cond_ind, :]

    @property
    def ires(self):
        return self.base_object.ires[self._cond_ind]

    @property
    def fwidth(self):
        return self.base_object.fwidth[self._cond_ind, :]

    def _get_bbox(self):
        """
        Get the bounding box for the cut region. Here we just use
        the bounding box for the source region.
        """
        return self.base_object._get_bbox()

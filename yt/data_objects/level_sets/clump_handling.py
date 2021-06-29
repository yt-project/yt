import uuid

import numpy as np

from yt.fields.derived_field import ValidateSpatial
from yt.frontends.ytdata.utilities import save_as_dataset
from yt.funcs import get_output_filename, mylog
from yt.utilities.tree_container import TreeContainer

from .clump_info_items import clump_info_registry
from .clump_validators import clump_validator_registry
from .contour_finder import identify_contours


def add_contour_field(ds, contour_key):
    def _contours(field, data):
        fd = data.get_field_parameter(f"contour_slices_{contour_key}")
        vals = data["index", "ones"] * -1
        if fd is None or fd == 0.0:
            return vals
        for sl, v in fd.get(data.id, []):
            vals[sl] = v
        return vals

    ds.add_field(
        ("index", f"contours_{contour_key}"),
        function=_contours,
        validators=[ValidateSpatial(0)],
        take_log=False,
        display_field=False,
        sampling_type="cell",
        units="",
    )


class Clump(TreeContainer):
    def __init__(
        self,
        data,
        field,
        parent=None,
        clump_info=None,
        validators=None,
        base=None,
        contour_key=None,
        contour_id=None,
    ):
        self.data = data
        self.field = field
        self.parent = parent
        self.quantities = data.quantities
        self.min_val = self.data[field].min()
        self.max_val = self.data[field].max()
        self.info = {}
        self.children = []

        # is this the parent clump?
        if base is None:
            base = self
            self.total_clumps = 0

        if clump_info is None:
            self.set_default_clump_info()
        else:
            self.clump_info = clump_info

        for ci in self.clump_info:
            ci(self)

        self.base = base
        self.clump_id = self.base.total_clumps
        self.base.total_clumps += 1
        self.contour_key = contour_key
        self.contour_id = contour_id

        if parent is not None:
            self.data.parent = self.parent.data

        if validators is None:
            validators = []
        self.validators = validators
        # Return value of validity function.
        self.valid = None

    _leaves = None

    @property
    def leaves(self):
        if self._leaves is not None:
            return self._leaves

        self._leaves = []
        for clump in self:
            if not clump.children:
                self._leaves.append(clump)
        return self._leaves

    def add_validator(self, validator, *args, **kwargs):
        """
        Add a validating function to determine whether the clump should
        be kept.
        """
        callback = clump_validator_registry.find(validator, *args, **kwargs)
        self.validators.append(callback)
        for child in self.children:
            child.add_validator(validator)

    def add_info_item(self, info_item, *args, **kwargs):
        "Adds an entry to clump_info list and tells children to do the same."

        callback = clump_info_registry.find(info_item, *args, **kwargs)
        callback(self)
        self.clump_info.append(callback)
        for child in self.children:
            child.add_info_item(info_item)

    def set_default_clump_info(self):
        "Defines default entries in the clump_info array."

        # add_info_item is recursive so this function does not need to be.
        self.clump_info = []

        self.add_info_item("total_cells")
        self.add_info_item("cell_mass")

        if any("jeans" in f for f in self.data.pf.field_list):
            self.add_info_item("mass_weighted_jeans_mass")
            self.add_info_item("volume_weighted_jeans_mass")

        self.add_info_item("max_grid_level")

        if any("number_density" in f for f in self.data.pf.field_list):
            self.add_info_item("min_number_density")
            self.add_info_item("max_number_density")

    def clear_clump_info(self):
        """
        Clears the clump_info array and passes the instruction to its
        children.
        """

        self.clump_info = []
        for child in self.children:
            child.clear_clump_info()

    def find_children(self, min_val, max_val=None):
        if self.children:
            mylog.info("Wiping out existing children clumps: %d.", len(self.children))
        self.children = []
        if max_val is None:
            max_val = self.max_val
        nj, cids = identify_contours(self.data, self.field, min_val, max_val)
        # Here, cids is the set of slices and values, keyed by the
        # parent_grid_id, that defines the contours.  So we can figure out all
        # the unique values of the contours by examining the list here.
        unique_contours = set()
        for sl_list in cids.values():
            for _sl, ff in sl_list:
                unique_contours.update(np.unique(ff))
        contour_key = uuid.uuid4().hex
        base_object = getattr(self.data, "base_object", self.data)
        add_contour_field(base_object.ds, contour_key)
        for cid in sorted(unique_contours):
            if cid == -1:
                continue
            new_clump = base_object.cut_region(
                [f"obj['contours_{contour_key}'] == {cid}"],
                {(f"contour_slices_{contour_key}"): cids},
            )
            if new_clump[("index", "ones")].size == 0:
                # This is to skip possibly duplicate clumps.
                # Using "ones" here will speed things up.
                continue
            self.children.append(
                Clump(
                    new_clump,
                    self.field,
                    parent=self,
                    validators=self.validators,
                    base=self.base,
                    clump_info=self.clump_info,
                    contour_key=contour_key,
                    contour_id=cid,
                )
            )

    def __iter__(self):
        yield self
        for child in self.children:
            yield from child

    def save_as_dataset(self, filename=None, fields=None):
        r"""Export clump tree to a reloadable yt dataset.
        This function will take a clump object and output a dataset
        containing the fields given in the ``fields`` list and all info
        items.  The resulting dataset can be reloaded as a yt dataset.

        Parameters
        ----------
        filename : str, optional
            The name of the file to be written.  If None, the name
            will be a combination of the original dataset and the clump
            index.
        fields : list of strings or tuples, optional
            If this is supplied, it is the list of fields to be saved to
            disk.

        Returns
        -------
        filename : str
            The name of the file that has been created.

        Examples
        --------

        >>> import yt
        >>> from yt.data_objects.level_sets.api import Clump, find_clumps
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> data_source = ds.disk(
        ...     [0.5, 0.5, 0.5], [0.0, 0.0, 1.0], (8, "kpc"), (1, "kpc")
        ... )
        >>> field = ("gas", "density")
        >>> step = 2.0
        >>> c_min = 10 ** np.floor(np.log10(data_source[field]).min())
        >>> c_max = 10 ** np.floor(np.log10(data_source[field]).max() + 1)
        >>> master_clump = Clump(data_source, field)
        >>> master_clump.add_info_item("center_of_mass")
        >>> master_clump.add_validator("min_cells", 20)
        >>> find_clumps(master_clump, c_min, c_max, step)
        >>> fn = master_clump.save_as_dataset(
        ...     fields=[("gas", "density"), ("all", "particle_mass")]
        ... )
        >>> new_ds = yt.load(fn)
        >>> print(ds.tree["clump", "cell_mass"])
        1296926163.91 Msun
        >>> print(ds.tree["grid", "density"])
        [  2.54398434e-26   2.46620353e-26   2.25120154e-26 ...,   1.12879234e-25
           1.59561490e-25   1.09824903e-24] g/cm**3
        >>> print(ds.tree["all", "particle_mass"])
        [  4.25472446e+38   4.25472446e+38   4.25472446e+38 ...,   2.04238266e+38
           2.04523901e+38   2.04770938e+38] g
        >>> print(ds.tree.children[0]["clump", "cell_mass"])
        909636495.312 Msun
        >>> print(ds.leaves[0]["clump", "cell_mass"])
        3756566.99809 Msun
        >>> print(ds.leaves[0]["grid", "density"])
        [  6.97820274e-24   6.58117370e-24   7.32046082e-24   6.76202430e-24
           7.41184837e-24   6.76981480e-24   6.94287213e-24   6.56149658e-24
           6.76584569e-24   6.94073710e-24   7.06713082e-24   7.22556526e-24
           7.08338898e-24   6.78684331e-24   7.40647040e-24   7.03050456e-24
           7.12438678e-24   6.56310217e-24   7.23201662e-24   7.17314333e-24] g/cm**3

        """

        ds = self.data.ds
        keyword = "%s_clump_%d" % (str(ds), self.clump_id)
        filename = get_output_filename(filename, keyword, ".h5")

        # collect clump info fields
        clump_info = {ci.name: [] for ci in self.base.clump_info}
        clump_info.update(
            {
                field: []
                for field in ["clump_id", "parent_id", "contour_key", "contour_id"]
            }
        )
        for clump in self:
            clump_info["clump_id"].append(clump.clump_id)
            if clump.parent is None:
                parent_id = -1
            else:
                parent_id = clump.parent.clump_id
            clump_info["parent_id"].append(parent_id)

            contour_key = clump.contour_key
            if contour_key is None:
                contour_key = -1
            clump_info["contour_key"].append(contour_key)
            contour_id = clump.contour_id
            if contour_id is None:
                contour_id = -1
            clump_info["contour_id"].append(contour_id)

            for ci in self.base.clump_info:
                clump_info[ci.name].append(clump.info[ci.name][1])
        for ci in clump_info:
            if hasattr(clump_info[ci][0], "units"):
                clump_info[ci] = ds.arr(clump_info[ci])
            else:
                clump_info[ci] = np.array(clump_info[ci])

        ftypes = {ci: "clump" for ci in clump_info}

        # collect data fields
        if fields is not None:
            contour_fields = [
                ("index", f"contours_{ckey}")
                for ckey in np.unique(clump_info["contour_key"])
                if str(ckey) != "-1"
            ]

            ptypes = []
            field_data = {}
            need_grid_positions = False
            for f in self.base.data._determine_fields(fields) + contour_fields:
                if ds.field_info[f].sampling_type == "particle":
                    if f[0] not in ptypes:
                        ptypes.append(f[0])
                    ftypes[f] = f[0]
                else:
                    need_grid_positions = True
                    if f[1] in ("x", "y", "z", "dx", "dy", "dz"):
                        # skip 'xyz' if a user passes that in because they
                        # will be added to ftypes below
                        continue
                    ftypes[f] = "grid"
                field_data[f] = self.base[f]

            if len(ptypes) > 0:
                for ax in "xyz":
                    for ptype in ptypes:
                        p_field = (ptype, f"particle_position_{ax}")
                        if p_field in ds.field_info and p_field not in field_data:
                            ftypes[p_field] = p_field[0]
                            field_data[p_field] = self.base[p_field]

                for clump in self:
                    if clump.contour_key is None:
                        continue
                    for ptype in ptypes:
                        cfield = (ptype, f"contours_{clump.contour_key}")
                        if cfield not in field_data:
                            field_data[cfield] = clump.data._part_ind(ptype).astype(
                                np.int64
                            )
                            ftypes[cfield] = ptype
                        field_data[cfield][
                            clump.data._part_ind(ptype)
                        ] = clump.contour_id

            if need_grid_positions:
                for ax in "xyz":
                    g_field = ("index", ax)
                    if g_field in ds.field_info and g_field not in field_data:
                        field_data[g_field] = self.base[g_field]
                        ftypes[g_field] = "grid"
                    g_field = ("index", "d" + ax)
                    if g_field in ds.field_info and g_field not in field_data:
                        ftypes[g_field] = "grid"
                        field_data[g_field] = self.base[g_field]

            if self.contour_key is not None:
                cfilters = {}
                for field in field_data:
                    if ftypes[field] == "grid":
                        ftype = "index"
                    else:
                        ftype = field[0]
                    cfield = (ftype, f"contours_{self.contour_key}")
                    if cfield not in cfilters:
                        cfilters[cfield] = field_data[cfield] == self.contour_id
                    field_data[field] = field_data[field][cfilters[cfield]]

            clump_info.update(field_data)
        extra_attrs = {"data_type": "yt_clump_tree", "container_type": "yt_clump_tree"}
        save_as_dataset(
            ds, filename, clump_info, field_types=ftypes, extra_attrs=extra_attrs
        )

        return filename

    def pass_down(self, operation):
        """
        Performs an operation on a clump with an exec and passes the
        instruction down to clump children.
        """

        # Call if callable, otherwise do an exec.
        if callable(operation):
            operation()
        else:
            exec(operation)

        for child in self.children:
            child.pass_down(operation)

    def _validate(self):
        "Apply all user specified validator functions."

        # Only call functions if not done already.
        if self.valid is not None:
            return self.valid

        self.valid = True
        for validator in self.validators:
            self.valid &= validator(self)
            if not self.valid:
                break

        return self.valid

    def __reduce__(self):
        raise RuntimeError(
            "Pickling Clump instances is not supported. Please use "
            "Clump.save_as_dataset instead"
        )

    def __getitem__(self, request):
        return self.data[request]


def find_clumps(clump, min_val, max_val, d_clump):
    mylog.info("Finding clumps: min: %e, max: %e, step: %f", min_val, max_val, d_clump)
    if min_val >= max_val:
        return
    clump.find_children(min_val, max_val=max_val)

    if len(clump.children) == 1:
        find_clumps(clump, min_val * d_clump, max_val, d_clump)

    elif len(clump.children) > 0:
        these_children = []
        mylog.info("Investigating %d children.", len(clump.children))
        for child in clump.children:
            find_clumps(child, min_val * d_clump, max_val, d_clump)
            if len(child.children) > 0:
                these_children.append(child)
            elif child._validate():
                these_children.append(child)
            else:
                mylog.info(
                    "Eliminating invalid, childless clump with %d cells.",
                    len(child.data[("index", "ones")]),
                )
        if len(these_children) > 1:
            mylog.info(
                "%d of %d children survived.", len(these_children), len(clump.children)
            )
            clump.children = these_children
        elif len(these_children) == 1:
            mylog.info(
                "%d of %d children survived, linking its children to parent.",
                len(these_children),
                len(clump.children),
            )
            clump.children = these_children[0].children
            for child in clump.children:
                child.parent = clump
                child.data.parent = clump.data
        else:
            mylog.info(
                "%d of %d children survived, erasing children.",
                len(these_children),
                len(clump.children),
            )
            clump.children = []

"""
Data structures for SPH frontends.




"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.data_objects.static_output import \
    Dataset

class ParticleDataset(Dataset):
    _unit_base = None
    over_refine_factor = 1
    filter_bbox = False

    def add_deposited_particle_field(self, deposit_field, method, kernel_name='cubic',
                                     weight_field='particle_mass'):
        """Add a new deposited particle field

        Creates a new deposited field based on the particle *deposit_field*.

        Parameters
        ----------

        deposit_field : tuple
           The field name tuple of the particle field the deposited field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        method : string
           This is the "method name" which will be looked up in the
           `particle_deposit` namespace as `methodname_deposit`.  Current
           methods include `simple_smooth`, `sum`, `std`, `cic`, `weighted_mean`,
           `mesh_id`, and `nearest`.
        kernel_name : string, default 'cubic'
           This is the name of the smoothing kernel to use. It is only used for
           the `simple_smooth` method and is otherwise ignored. Current
           supported kernel names include `cubic`, `quartic`, `quintic`,
           `wendland2`, `wendland4`, and `wendland6`.
        weight_field : string, default 'particle_mass'
           Weighting field name for deposition method `weighted_mean`.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        self.index
        if isinstance(deposit_field, tuple):
            ptype, deposit_field = deposit_field[0], deposit_field[1]
        else:
            raise RuntimeError

        units = self.field_info[ptype, deposit_field].units
        take_log = self.field_info[ptype, deposit_field].take_log
        name_map = {"sum": "sum", "std":"std", "cic": "cic", "weighted_mean": "avg",
                    "nearest": "nn", "simple_smooth": "ss", "count": "count"}
        field_name = "%s_" + name_map[method] + "_%s"
        field_name = field_name % (ptype, deposit_field.replace('particle_', ''))

        if method == "count":
            field_name = "%s_count" % ptype
            if ("deposit", field_name) in self.field_info:
                mylog.warning("The deposited field %s already exists" % field_name)
                return ("deposit", field_name)
            else:
                units = "dimensionless"
                take_log = False

        def _deposit_field(field, data):
            """
            Create a grid field for particle quantities using given method.
            """
            pos = data[ptype, "particle_position"]
            if method == 'weighted_mean':
                d = data.ds.arr(data.deposit(pos, [data[ptype, deposit_field],
                                                   data[ptype, weight_field]],
                                             method=method, kernel_name=kernel_name),
                                             input_units=units)
                d[np.isnan(d)] = 0.0
            else:
                d = data.ds.arr(data.deposit(pos, [data[ptype, deposit_field]],
                                             method=method, kernel_name=kernel_name),
                                             input_units=units)
            return d

        self.add_field(
            ("deposit", field_name),
            function=_deposit_field,
            units=units,
            take_log=take_log,
            validators=[ValidateSpatial()])
        return ("deposit", field_name)

    def add_smoothed_particle_field(self, smooth_field, method="volume_weighted",
                                    nneighbors=64, kernel_name="cubic"):
        """Add a new smoothed particle field

        Creates a new smoothed field based on the particle *smooth_field*.

        Parameters
        ----------

        smooth_field : tuple
           The field name tuple of the particle field the smoothed field will
           be created from.  This must be a field name tuple so yt can
           appropriately infer the correct particle type.
        method : string, default 'volume_weighted'
           The particle smoothing method to use. Can only be 'volume_weighted'
           for now.
        nneighbors : int, default 64
            The number of neighbors to examine during the process.
        kernel_name : string, default 'cubic'
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        self.index
        if isinstance(smooth_field, tuple):
            ptype, smooth_field = smooth_field[0], smooth_field[1]
        else:
            raise RuntimeError("smooth_field must be a tuple, received %s" %
                               smooth_field)
        if method != "volume_weighted":
            raise NotImplementedError("method must be 'volume_weighted'")

        coord_name = "particle_position"
        mass_name = "particle_mass"
        smoothing_length_name = "smoothing_length"
        if (ptype, smoothing_length_name) not in self.derived_field_list:
            raise ValueError("%s not in derived_field_list" %
                             ((ptype, smoothing_length_name),))
        density_name = "density"
        registry = self.field_info

        return add_volume_weighted_smoothed_field(ptype, coord_name, mass_name,
                   smoothing_length_name, density_name, smooth_field, registry,
                   nneighbors=nneighbors, kernel_name=kernel_name)[0]

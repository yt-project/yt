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
from yt.fields.particle_fields import \
    add_volume_weighted_smoothed_field


class ParticleDataset(Dataset):
    _unit_base = None
    default_kernel = "cubic"
    filter_bbox = False
    over_refine_factor = 1

    def add_smoothed_particle_field(self, smooth_field,
                                    method="volume_weighted", nneighbors=64,
                                    kernel_name=None):
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
        kernel_name : string or None, default None
            This is the name of the smoothing kernel to use. Current supported
            kernel names include `cubic`, `quartic`, `quintic`, `wendland2`,
            `wendland4`, and `wendland6`. If left as None,
            :attr:`~yt.frontends.sph.data_structures.ParticleDataset.default_kernel`
            will be used.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        # The magical step
        self.index

        # Parse arguments
        if isinstance(smooth_field, tuple):
            ptype, smooth_field = smooth_field[0], smooth_field[1]
        else:
            raise RuntimeError("smooth_field must be a tuple, received %s" %
                               smooth_field)
        if method != "volume_weighted":
            raise NotImplementedError("method must be 'volume_weighted'")
        if kernel_name is None:
            kernel_name = self.default_kernel

        # Prepare field names and registry to be used later
        coord_name = "particle_position"
        mass_name = "particle_mass"
        smoothing_length_name = "smoothing_length"
        if (ptype, smoothing_length_name) not in self.derived_field_list:
            raise ValueError("%s not in derived_field_list" %
                             ((ptype, smoothing_length_name),))
        density_name = "density"
        registry = self.field_info

        # Do the actual work
        return add_volume_weighted_smoothed_field(ptype, coord_name, mass_name,
                   smoothing_length_name, density_name, smooth_field, registry,
                   nneighbors=nneighbors, kernel_name=kernel_name)[0]

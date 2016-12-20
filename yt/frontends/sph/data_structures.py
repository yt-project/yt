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
    ParticleDataset


class SPHDataset(ParticleDataset):
    default_kernel_name = "cubic"

    def __init__(self, filename, dataset_type=None, file_style=None,
                 units_override=None, unit_system="cgs",
                 n_ref=64, over_refine_factor=1,
                 kernel_name=None):
        if kernel_name is None:
            self.kernel_name = self.default_kernel_name
        else:
            self.kernel_name = kernel_name
        super(SPHDataset, self).__init__(
            filename, dataset_type=dataset_type, file_style=file_style,
            units_override=units_override, unit_system=unit_system,
            n_ref=n_ref, over_refine_factor=over_refine_factor)

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
            :attr:`~yt.frontends.sph.data_structures.SPHDataset.kernel_name`
            will be used.

        Returns
        -------

        The field name tuple for the newly created field.
        """
        if kernel_name is None:
            kernel_name = self.kernel_name
        return super(SPHDataset, self).add_smoothed_particle_field(
            smooth_field=smooth_field, method=method, nneighbors=nneighbors,
            kernel_name=kernel_name
        )

"""
Gadget-specfic fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.sph.fields import SPHFieldInfo
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
from yt.utilities.physical_constants import mp, kb

class GadgetFieldInfo(SPHFieldInfo):

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super(GadgetFieldInfo, self).setup_particle_fields(
            ptype, *args, **kwargs)

        # setup some special fields that only make sense for SPH particles
        if ptype in ("PartType0", "Gas"):
            self.setup_gas_particle_fields(ptype)

    def setup_gas_particle_fields(self, ptype):
        if (ptype, "ElectronAbundance") in self.ds.field_list:
            def _temperature(field, data):
                # Assume cosmic abundances
                x_H = 0.76
                gamma = 5.0/3.0
                a_e = data[ptype, 'ElectronAbundance']
                mu = 4.0 / (3.0 * x_H + 1.0 + 4.0 * x_H * a_e)
                ret = data[ptype, "InternalEnergy"]*(gamma-1)*mu*mp/kb
                return ret.in_units('K')
        else:
            def _temperature(field, data):
                # Assume cosmic abundances
                x_H = 0.76
                gamma = 5.0/3.0
                # Assume zero ionization
                mu = 4.0 / (3.0 * x_H + 1.0)
                ret = data[ptype, "InternalEnergy"]*(gamma-1)*mu*mp/kb
                return ret.in_units('K')

        self.add_field(
            (ptype, "Temperature"),
            function=_temperature,
            particle_type=True,
            units="K")

        # For now, we hardcode num_neighbors.  We should make this configurable
        # in the future.
        num_neighbors = 64
        fn = add_volume_weighted_smoothed_field(
            ptype, "particle_position", "particle_mass", "smoothing_length",
            "density", "Temperature", self, num_neighbors)

        # Alias ("gas", "temperature") to the new smoothed Temperature field
        self.alias(("gas", "temperature"), fn[0])

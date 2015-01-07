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
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.physical_constants import mp, kb

class GadgetFieldInfo(SPHFieldInfo):

    def setup_particle_fields(self, ptype, *args, **kwargs):
        super(GadgetFieldInfo, self).setup_particle_fields(
            ptype, *args, **kwargs)

        # Only add the Temperature field for SPH particles
        if ptype not in ("PartType0", "Gas"):
            return

        def _temperature(field, data):
            # Assume cosmic abundances
            x_H = 0.76
            gamma = 5.0/3.0
            try:
                a_e = data['ElectronAbundance']
            except YTFieldNotFound:
                # Assume zero ionization
                a_e = 0.0
            mu = 4.0 / (3.0 * x_H + 1.0 + 4.0 * x_H * a_e)
            ret = data[ptype, "InternalEnergy"]*(gamma-1)*mu*mp/kb
            return ret.in_units('K')

        self.add_field(
            (ptype, "Temperature"),
            function=_temperature,
            particle_type=True,
            units="K")

        num_neighbors = 64
        fn = add_volume_weighted_smoothed_field(
            ptype, "particle_position", "particle_mass", "smoothing_length",
            "density", "Temperature", self, num_neighbors)

        # Alias ("gas", "temperature") to the new smoothed Temperature field
        self.alias(("gas", "temperature"), fn[0])

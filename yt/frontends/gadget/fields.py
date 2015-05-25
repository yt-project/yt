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

        # setup some special fields that only make sense for SPH particles
        if (ptype, "FourMetalFractions") in self.ds.field_list:
            metal_fraction_names = ['C', 'O', 'Si', 'Fe']
            for i in range(len(metal_fraction_names)):
                def _Fraction_wrap(i):
                    def _Fraction(field,data):
                        return data[(ptype, 'FourMetalFractions')][:,i]
                    return _Fraction
#                def _Deposit_wrap(i):
#                    def _deposit(field, data):
#                        pos = data[ptype, coord_name]
#                        mass = data[ptype, mass_name]
#                        pos.convert_to_units("code_length")
#                        mass.convert_to_units("code_mass")
#                        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
#                        d = data.ds.arr(d, "code_mass")
#                        d /= data["index", "cell_volume"]
#                        return d 
#                    return _deposit

                self.add_field( (ptype, metal_fraction_names[i]+"_fraction"),
                                function=_Fraction_wrap(i), 
                                particle_type=True,
                                units="")

        super(GadgetFieldInfo, self).setup_particle_fields(
            ptype, *args, **kwargs)

#                self.add_field( ("deposit", "%s_%s" % (ptype, metal_fraction_names[i]+"Mass")),
#                                function=_Deposit_wrap(i), 
#                                units="")

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

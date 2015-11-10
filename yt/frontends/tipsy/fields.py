"""
Tipsy fields




"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.sph.fields import SPHFieldInfo
from yt.fields.particle_fields import add_nearest_neighbor_field

class TipsyFieldInfo(SPHFieldInfo):
    aux_particle_fields = {
        'uDotFB':("uDotFB", ("code_mass * code_velocity**2", [""], None)),
        'uDotAV':("uDotAV", ("code_mass * code_velocity**2", [""], None)),
        'uDotPdV':("uDotPdV", ("code_mass * code_velocity**2", [""], None)),
        'uDotHydro':("uDotHydro", ("code_mass * code_velocity**2", [""], None)),
        'uDotDiff':("uDotDiff", ("code_mass * code_velocity**2", [""], None)),
        'uDot':("uDot", ("code_mass * code_velocity**2", [""], None)),
        'coolontime':("coolontime", ("code_time", [""], None)),
        'timeform':("timeform", ("code_time", [""], None)),
        'massform':("massform", ("code_mass", [""], None)),
        'HI':("HI", ("dimensionless", ["H_fraction"], None)),
        'HII':("HII", ("dimensionless", ["H_p1_fraction"], None)),
        'HeI':("HeI", ("dimensionless", ["He_fraction"], None)),
        'HeII':("HeII", ("dimensionless", ["He_p2_fraction"], None)),
        'OxMassFrac':("OxMassFrac", ("dimensionless", ["O_fraction"], None)),
        'FeMassFrac':("FeMassFrac", ("dimensionless", ["Fe_fraction"], None)),
        'c':("c", ("code_velocity", [""], None)),
        'acc':("acc", ("code_velocity / code_time", [""], None)),
        'accg':("accg", ("code_velocity / code_time", [""], None)),
        'smoothlength':('smoothlength', ("code_length", ["smoothing_length"], None))}

    def __init__(self, ds, field_list, slice_info = None):
        for field in field_list:
            if field[1] in self.aux_particle_fields.keys() and \
                self.aux_particle_fields[field[1]] not in self.known_particle_fields:
                self.known_particle_fields += (self.aux_particle_fields[field[1]],)
        super(TipsyFieldInfo,self).__init__(ds, field_list, slice_info)

    def setup_particle_fields(self, ptype, *args, **kwargs):

        # setup some special fields that only make sense for SPH particles

        if ptype in ("PartType0", "Gas"):
            self.setup_gas_particle_fields(ptype)

        super(TipsyFieldInfo, self).setup_particle_fields(
            ptype, *args, **kwargs)


    def setup_gas_particle_fields(self, ptype):

        num_neighbors = 65
        fn, = add_nearest_neighbor_field(ptype, "particle_position", self, num_neighbors)
        def _func():
            def _smoothing_length(field, data):
                # For now, we hardcode num_neighbors.  We should make this configurable
                # in the future.
                rv = data[ptype, 'nearest_neighbor_distance_%d' % num_neighbors]
                #np.maximum(rv, 0.5*data[ptype, "Epsilon"], rv)
                return rv
            return _smoothing_length

        self.add_field(
            (ptype, "smoothing_length"),
            function=_func(),
            particle_type=True,
            units="code_length")

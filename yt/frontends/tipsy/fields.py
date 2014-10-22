"""
Tipsy fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.sph.fields import \
    SPHFieldInfo

class TipsyFieldInfo(SPHFieldInfo):
    aux_particle_fields = {
        'uDotFB':("uDotFB", ("code_mass * code_velocity**2", ["uDotFB"], None)),
        'uDotAV':("uDotAV", ("code_mass * code_velocity**2", ["uDotAV"], None)),
        'uDotPdV':("uDotPdV", ("code_mass * code_velocity**2", ["uDotPdV"], None)),
        'uDotHydro':("uDotHydro", ("code_mass * code_velocity**2", ["uDotHydro"], None)),
        'uDotDiff':("uDotDiff", ("code_mass * code_velocity**2", ["uDotDiff"], None)),
        'uDot':("uDot", ("code_mass * code_velocity**2", ["uDot"], None)),
        'coolontime':("coolontime", ("code_time", ["coolontime"], None)),
        'timeform':("timeform", ("code_time", ["timeform"], None)),
        'massform':("massform", ("code_mass", ["massform"], None)),
        'HI':("HI", ("dimensionless", ["HI"], None)),
        'HII':("HII", ("dimensionless", ["HII"], None)),
        'HeI':("HeI", ("dimensionless", ["HeI"], None)),
        'HeII':("HeII", ("dimensionless", ["HeII"], None)),
        'OxMassFrac':("OxMassFrac", ("dimensionless", ["OxMassFrac"], None)),
        'FeMassFrac':("FeMassFrac", ("dimensionless", ["FeMassFrac"], None)),
        'c':("c", ("code_velocity", ["c"], None)),
        'acc':("acc", ("code_velocity / code_time", ["acc"], None)),
        'accg':("accg", ("code_velocity / code_time", ["accg"], None))}
    
    def __init__(self, ds, field_list, slice_info = None):
        for field in field_list:
            if field[1] in self.aux_particle_fields.keys() and \
                self.aux_particle_fields[field[1]] not in self.known_particle_fields:
                self.known_particle_fields += (self.aux_particle_fields[field[1]],)
        super(TipsyFieldInfo,self).__init__(ds, field_list, slice_info)

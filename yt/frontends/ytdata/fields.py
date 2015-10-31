"""
YTData-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer

m_units = "g"
p_units = "cm"
v_units = "cm / s"
r_units = "cm"

class YTDataContainerFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("x", (p_units, ["particle_position_x"], None)),
        ("y", (p_units, ["particle_position_y"], None)),
        ("z", (p_units, ["particle_position_z"], None)),
        ("velocity_x", (v_units, ["particle_velocity_x"], None)),
        ("velocity_y", (v_units, ["particle_velocity_y"], None)),
        ("velocity_z", (v_units, ["particle_velocity_z"], None)),
    )

    # these are extra fields to be created for the "all" particle type
    extra_union_fields = (
        (p_units, "particle_position_x"),
        (p_units, "particle_position_y"),
        (p_units, "particle_position_z"),
        (v_units, "particle_velocity_x"),
        (v_units, "particle_velocity_y"),
        (v_units, "particle_velocity_z"),
    )

    def __init__(self, ds, field_list):
        super(YTDataContainerFieldInfo, self).__init__(ds, field_list)
        self.add_fake_grid_fields()

    def add_fake_grid_fields(self):
        """
        Add cell volume and mass fields that use the dx, dy, and dz
        fields that come with the dataset instead of the index fields
        which correspond to the oct tree.  We need to do this for now
        since we're treating the grid data like particles until we
        implement exporting AMR hierarchies.
        """

        if ("grid", "cell_volume") not in self.field_list:
            def _cell_volume(field, data):
                return data["grid", "dx"] * \
                  data["grid", "dy"] * \
                  data["grid", "dz"]
            self.add_field(("grid", "cell_volume"), function=_cell_volume,
                           units="cm**3", particle_type=True)

class YTGridFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
    )

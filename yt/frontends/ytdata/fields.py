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

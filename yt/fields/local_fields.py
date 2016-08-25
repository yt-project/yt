"""
This is a container for storing local fields defined on each load of yt.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.logger import \
    ytLogger as mylog

from .field_plugin_registry import \
    register_field_plugin

from .field_info_container import \
    FieldInfoContainer

class LocalFieldInfoContainer(FieldInfoContainer):
    def add_field(self, name, function=None, **kwargs):
        if not isinstance(name, tuple):
            if kwargs.setdefault('particle_type', False):
                name = ('all', name)
            else:
                name = ('gas', name)
        override = kwargs.get("force_override", False)
        # Handle the case where the field has already been added.
        if not override and name in self:
            mylog.warning("Field %s already exists. To override use " +
                          "force_override=True.", name)
        return super(LocalFieldInfoContainer,
                     self).add_field(name, function, **kwargs)

# Empty FieldInfoContainer
local_fields = LocalFieldInfoContainer(None, [], None)

add_field = derived_field = local_fields.add_field

@register_field_plugin
def setup_local_fields(registry, ftype = "gas", slice_info = None):
    # This is easy.  We just update with the contents of the local_fields field
    # info container, and since they are not mutable in any real way, we are
    # fine.
    # Note that we actually don't care about the ftype here.
    for f in local_fields:
        registry._show_field_errors.append(f)
    registry.update(local_fields)

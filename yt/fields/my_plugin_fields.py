"""
This is a container for storing fields defined in the my_plugins.py file.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from .field_plugin_registry import \
    register_field_plugin

from .field_info_container import \
    FieldInfoContainer

# Empty FieldInfoContainer
my_plugins_fields = FieldInfoContainer(None, [], None)

@register_field_plugin
def setup_my_plugins_fields(registry, ftype="gas", slice_info=None):
    # fields end up inside this container when added via add_field in
    # my_plugins.py. See yt.funcs.enable_plugins to see how this is set up.
    registry.update(my_plugins_fields)

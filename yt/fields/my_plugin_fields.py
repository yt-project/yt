from typing import Any, Optional

from yt.frontends.ramses.fields import RAMSESFieldInfo

from .field_plugin_registry import register_field_plugin
from .local_fields import LocalFieldInfoContainer

# Empty FieldInfoContainer
my_plugins_fields = LocalFieldInfoContainer(None, [], None)


@register_field_plugin
def setup_my_plugins_fields(
    registry: RAMSESFieldInfo, ftype: str = "gas", slice_info: Optional[Any] = None
) -> None:
    # fields end up inside this container when added via add_field in
    # my_plugins.py. See yt.funcs.enable_plugins to see how this is set up.
    registry.update(my_plugins_fields)

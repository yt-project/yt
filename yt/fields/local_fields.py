from yt.funcs import is_sequence
from yt.utilities.logger import ytLogger as mylog

from .field_info_container import FieldInfoContainer
from .field_plugin_registry import register_field_plugin


class LocalFieldInfoContainer(FieldInfoContainer):
    def add_field(self, name, function, sampling_type, **kwargs):

        sampling_type = self._sanitize_sampling_type(
            sampling_type, kwargs.get("particle_type")
        )

        if isinstance(name, str) or not is_sequence(name):
            if sampling_type == "particle":
                ftype = "all"
            else:
                ftype = "gas"
            name = (ftype, name)

        override = kwargs.get("force_override", False)
        # Handle the case where the field has already been added.
        if not override and name in self:
            mylog.warning(
                "Field %s already exists. To override use `force_override=True`.",
                name,
            )

        return super().add_field(name, function, sampling_type, **kwargs)


# Empty FieldInfoContainer
local_fields = LocalFieldInfoContainer(None, [], None)

# we define two handles, pointing to the same function but documented differently
# yt.add_field() is meant to be used directly, while yt.derived_field is documented
# as a decorator.
add_field = local_fields.add_field


def derived_field(name, sampling_type, **kwargs):
    return add_field(name=name, function=None, sampling_type=sampling_type, **kwargs)


@register_field_plugin
def setup_local_fields(registry, ftype="gas", slice_info=None):
    # This is easy.  We just update with the contents of the local_fields field
    # info container, and since they are not mutable in any real way, we are
    # fine.
    # Note that we actually don't care about the ftype here.
    for f in local_fields:
        registry._show_field_errors.append(f)
    registry.update(local_fields)

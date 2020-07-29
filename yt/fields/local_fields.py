import warnings

from yt.utilities.logger import ytLogger as mylog

from .field_info_container import FieldInfoContainer
from .field_plugin_registry import register_field_plugin


class LocalFieldInfoContainer(FieldInfoContainer):
    def add_field(self, name, function=None, sampling_type=None, **kwargs):
        if not isinstance(name, tuple):
            if kwargs.setdefault("particle_type", False):
                name = ("all", name)
            else:
                name = ("gas", name)
        override = kwargs.get("force_override", False)
        # Handle the case where the field has already been added.
        if not override and name in self:
            mylog.error(
                "Field %s already exists. To override use " + "force_override=True.",
                name,
            )
        if kwargs.setdefault("particle_type", False):
            if sampling_type is not None and sampling_type != "particle":
                raise RuntimeError(
                    "Clashing definition of 'sampling_type' and "
                    "'particle_type'. Note that 'particle_type' is "
                    "deprecated. Please just use 'sampling_type'."
                )
            else:
                sampling_type = "particle"
        if sampling_type is None:
            warnings.warn(
                "Because 'sampling_type' is not specified, yt will "
                "assume a 'cell' sampling_type for the %s field" % (name,),
                stacklevel=3,
            )
            sampling_type = "cell"
        return super(LocalFieldInfoContainer, self).add_field(
            name, sampling_type, function, **kwargs
        )


# Empty FieldInfoContainer
local_fields = LocalFieldInfoContainer(None, [], None)

add_field = derived_field = local_fields.add_field


@register_field_plugin
def setup_local_fields(registry, ftype="gas", slice_info=None):
    # This is easy.  We just update with the contents of the local_fields field
    # info container, and since they are not mutable in any real way, we are
    # fine.
    # Note that we actually don't care about the ftype here.
    for f in local_fields:
        registry._show_field_errors.append(f)
    registry.update(local_fields)

from functools import partial
from typing import Any, Callable, TypeVar

from yt.funcs import is_sequence
from yt.utilities.logger import ytLogger as mylog

from .field_info_container import FieldInfoContainer
from .field_plugin_registry import register_field_plugin

# workaround mypy not being confortable around decorator preserving signatures
# adapted from
# https://github.com/python/mypy/issues/1551#issuecomment-253978622
TFun = TypeVar("TFun", bound=Callable[..., Any])


class LocalFieldInfoContainer(FieldInfoContainer):
    def add_field(
        self, name, function, sampling_type, *, force_override=False, **kwargs
    ):
        from yt.fields.field_functions import validate_field_function

        validate_field_function(function)
        if isinstance(name, str) or not is_sequence(name):
            # the base method only accepts proper tuple field keys
            # and is only used internally, while this method is exposed to users
            # and is documented as usable with single strings as name
            if sampling_type == "particle":
                ftype = "all"
            else:
                ftype = "gas"
            name = (ftype, name)

        # Handle the case where the field has already been added.
        if not force_override and name in self:
            mylog.warning(
                "Field %s already exists. To override use `force_override=True`.",
                name,
            )

        return super().add_field(
            name, function, sampling_type, force_override=force_override, **kwargs
        )


# Empty FieldInfoContainer
local_fields = LocalFieldInfoContainer(None, [], None)

# we define two handles, essentially pointing to the same function but documented differently
# yt.add_field() is meant to be used directly, while yt.derived_field is documented
# as a decorator.
add_field = local_fields.add_field


class derived_field:
    # implement a decorator accepting keyword arguments to be passed down to add_field

    def __init__(self, **kwargs) -> None:
        self._kwargs = kwargs

    def __call__(self, f: Callable) -> Callable:
        partial(local_fields.add_field, function=f)(**self._kwargs)
        return f


@register_field_plugin
def setup_local_fields(registry, ftype="gas", slice_info=None):
    # This is easy.  We just update with the contents of the local_fields field
    # info container, and since they are not mutable in any real way, we are
    # fine.
    # Note that we actually don't care about the ftype here.
    for f in local_fields:
        registry._show_field_errors.append(f)
    registry.update(local_fields)

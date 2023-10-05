from collections.abc import Callable

FunctionName = str
FieldPluginMap = dict[FunctionName, Callable]
field_plugins: FieldPluginMap = {}


def register_field_plugin(func: Callable) -> Callable:
    name = func.__name__
    if name.startswith("setup_"):
        name = name[len("setup_") :]
    if name.endswith("_fields"):
        name = name[: -len("_fields")]
    field_plugins[name] = func
    # And, we return it, too
    return func

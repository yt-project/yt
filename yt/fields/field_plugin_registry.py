field_plugins = {}


def register_field_plugin(func):
    name = func.__name__
    if name.startswith("setup_"):
        name = name[len("setup_") :]
    if name.endswith("_fields"):
        name = name[: -len("_fields")]
    field_plugins[name] = func
    # And, we return it, too
    return func

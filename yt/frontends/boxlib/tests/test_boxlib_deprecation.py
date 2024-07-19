from importlib import import_module, reload

from yt._maintenance.deprecation import warnings


def test_imports():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        for index, mname in enumerate(["data_structures", "fields", "io"]):
            mod_name = import_module("yt.frontends.boxlib." + mname)
            if len(w) != index + 1:
                reload(mod_name)

        assert len(w) == 3 and all(
            [
                issubclass(w[0].category, DeprecationWarning),
                issubclass(w[1].category, DeprecationWarning),
                issubclass(w[2].category, DeprecationWarning),
            ]
        )

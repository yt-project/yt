from yt._maintenance.deprecation import warnings


def test_imports():
    with warnings.catch_warnings(record=True) as w:
        from yt.frontends.boxlib import (
            data_structures,
            fields,
            io,
        )

        assert len(w) == 3 and all(
            [
                issubclass(w[0].category, DeprecationWarning),
                issubclass(w[1].category, DeprecationWarning),
                issubclass(w[2].category, DeprecationWarning),
            ]
        )
    del (data_structures, io, fields)
    warnings.resetwarnings()

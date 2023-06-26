__all__ = [
    "adaptahop",
    "ahf",
    "amrvac",
    "art",
    "arepo",
    "artio",
    "athena",
    "athena_pp",
    "boxlib",
    "cf_radial",
    "chimera",
    "chombo",
    "cholla",
    "eagle",
    "enzo_e",
    "enzo",
    "exodus_ii",
    "fits",
    "flash",
    "gadget",
    "gadget_fof",
    "gamer",
    "gdf",
    "gizmo",
    "halo_catalog",
    "http_stream",
    "moab",
    "nc4_cm1",
    "open_pmd",
    "owls",
    "owls_subfind",
    "ramses",
    "rockstar",
    "sdf",
    "stream",
    "swift",
    "tipsy",
    "ytdata",
]

from functools import lru_cache


@lru_cache(maxsize=None)
def __getattr__(value):
    import importlib

    if value == "_all":
        for _ in __all__:
            __getattr__(_)
        return

    if value not in __all__:
        raise AttributeError

    return importlib.import_module(f"yt.frontends.{value}.api")

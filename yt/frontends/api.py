import glob
import importlib
import os
import sys
import time
import types

_frontends = [
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


class _frontend_container:
    def __init__(self):
        for frontend in _frontends:
            _mod = f"yt.frontends.{frontend}.api"
            setattr(self, frontend, importlib.import_module(_mod))
        setattr(self, "api", importlib.import_module("yt.frontends.api"))
        setattr(self, "__name__", "yt.frontends.api")

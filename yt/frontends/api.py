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
    "chombo",
    "eagle",
    "enzo_p",
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
            _mod = "yt.frontends.%s.api" % frontend
            setattr(self, frontend, importlib.import_module(_mod))
        setattr(self, "api", importlib.import_module("yt.frontends.api"))
        setattr(self, "__name__", "yt.frontends.api")

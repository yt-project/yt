from yt.units.unit_registry import default_unit_registry
from unyt.unit_systems import add_constants

add_constants(globals(), registry=default_unit_registry)

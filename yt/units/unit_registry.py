from unyt.dimensions import dimensionless
from unyt.unit_registry import *

default_unit_registry = UnitRegistry(unit_system="cgs") # type: ignore

default_unit_registry.add("h", 1.0, dimensionless, tex_repr=r"h")

from unyt.unit_registry import *
from unyt.dimensions import dimensionless

default_unit_registry = UnitRegistry(unit_system='cgs')

default_unit_registry.add('h', 1.0, dimensionless, tex_repr=r"h")

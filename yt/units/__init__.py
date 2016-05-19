from yt.units import unit_symbols
from yt.utilities import physical_constants

from yt.units.yt_array import YTQuantity


# function to only import quantities into this namespace
# we go through the trouble of doing this instead of "import *"
# to avoid including extraneous variables (e.g. floating point
# constants used to *construct* a physical constant) in this namespace
def import_quantities(module, global_namespace):
    for key, value in module.__dict__.items():
        if isinstance(value, YTQuantity):
            global_namespace[key] = value

import_quantities(unit_symbols, globals())
import_quantities(physical_constants, globals())

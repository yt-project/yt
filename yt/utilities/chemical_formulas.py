"""
Very basic chemical formula parser.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import re
from .periodic_table import periodic_table

class ChemicalFormula:
    def __init__(self, formula_string):
        # See YTEP-0003 for information on the format.
        self.formula_string = formula_string
        self.elements = []
        if "_" in self.formula_string:
            molecule, ionization = self.formula_string.split("_")
            if ionization[0] == "p":
                charge = int(ionization[1:])
            elif ionization[0] == "m":
                charge = -int(ionization[1:])
            else:
                raise NotImplementedError
        elif self.formula_string.startswith('El'):
            molecule = self.formula_string
            charge = -1
        else:
            molecule = self.formula_string
            charge = 0
        self.charge = charge
        for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', molecule):
            if count == '': count = 1
            self.elements.append((periodic_table[element], int(count)))
        self.weight = sum(n * e.weight for e, n in self.elements)

    def __repr__(self):
        return self.formula_string

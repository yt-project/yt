import re

from .periodic_table import periodic_table
from .physical_ratios import _primordial_mass_fraction


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
        elif self.formula_string.startswith("El"):
            molecule = self.formula_string
            charge = -1
        else:
            molecule = self.formula_string
            charge = 0
        self.charge = charge
        for element, count in re.findall(r"([A-Z][a-z]*)(\d*)", molecule):
            if count == "":
                count = 1
            self.elements.append((periodic_table[element], int(count)))
        self.weight = sum(n * e.weight for e, n in self.elements)

    def __repr__(self):
        return self.formula_string


def compute_mu(ion_state):
    if ion_state == "ionized" or ion_state is None:
        # full ionization
        n_H = 2.0  # fully ionized hydrogen gives two particles
        n_He = 3.0  # fully ionized helium gives three particles
    elif ion_state == "neutral":
        n_H = 1.0  # neutral hydrogen gives one particle
        n_He = 1.0  # neutral helium gives one particle
    muinv = n_H * _primordial_mass_fraction["H"] / ChemicalFormula("H").weight
    muinv += n_He * _primordial_mass_fraction["He"] / ChemicalFormula("He").weight
    return 1.0 / muinv

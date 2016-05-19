from yt.testing import assert_equal
from yt.utilities.chemical_formulas import ChemicalFormula
from yt.utilities.periodic_table import periodic_table

_molecules = (

    ("H2O_p1", (("H", 2), ("O", 1)), 1),
    ("H2O_m1", (("H", 2), ("O", 1)), -1),
    ("H2O", (("H", 2), ("O", 1)), 0),
    ("H2SO4", (("H", 2), ("S", 1), ("O", 4)), 0),
    # Now a harder one
    ("UuoMtUuq3", (("Uuo", 1), ("Mt", 1), ("Uuq", 3)), 0)
)

def test_formulas():
    for formula, components, charge in _molecules:
        f = ChemicalFormula(formula)
        w = sum( n * periodic_table[e].weight for e, n in components)
        yield assert_equal, f.charge, charge
        yield assert_equal, f.weight, w
        for (n, c1), (e, c2) in zip(components, f.elements):
            yield assert_equal, n, e.symbol
            yield assert_equal, c1, c2

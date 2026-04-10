from numpy.testing import assert_allclose, assert_equal
import pytest

from yt.utilities.chemical_formulas import ChemicalFormula, compute_mu


def test_electron_formula_charge_weight_and_repr():
    electron = ChemicalFormula("El")
    ((element, count),) = electron.elements

    assert_equal(electron.charge, -1)
    assert_allclose(electron.weight, 0.00054858)
    assert_equal(repr(electron), "El")
    assert_equal(element.symbol, "El")
    assert_equal(count, 1)


def test_invalid_ionization_suffix_raises_not_implemented():
    with pytest.raises(NotImplementedError):
        ChemicalFormula("H2O_x1")


def test_compute_mu_unsupported_state_raises_unbound_local_error():
    with pytest.raises(UnboundLocalError):
        compute_mu("molecular")

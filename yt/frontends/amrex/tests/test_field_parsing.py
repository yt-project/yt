import pytest

from yt.frontends.amrex.fields import Substance


@pytest.mark.parametrize(
    "data, expected",
    [
        pytest.param("X(He5)", [("He", 5)], id="isotope_1"),
        pytest.param("X(C12)", [("C", 12)], id="isotope_2"),
        pytest.param("X(A1B2C3)", [("A", 1), ("B", 2), ("C", 3)], id="molecule_1"),
        pytest.param("X(C12H24)", [("C", 12), ("H", 24)], id="molecule_2"),
        pytest.param("X(H2O)", [("H", 2), ("O", 1)], id="molecule_3"),
        pytest.param("X(ash)", [("ash", 0)], id="descriptive_name"),
    ],
)
def test_Substance_spec(data, expected):
    assert Substance(data)._spec == expected


@pytest.mark.parametrize(
    "data, expected_type",
    [
        pytest.param("X(He5)", "isotope", id="isotope_1"),
        pytest.param("X(C12)", "isotope", id="isotope_2"),
        pytest.param("X(A1B2C3)", "molecule", id="molecule_1"),
        pytest.param("X(C12H24)", "molecule", id="molecule_2"),
        pytest.param("X(H2O)", "molecule", id="molecule_3"),
        pytest.param("X(ash)", "descriptive_name", id="descriptive_name"),
    ],
)
def test_Substance_type(data, expected_type):
    sub = Substance(data)
    assert getattr(sub, f"is_{expected_type}")()


@pytest.mark.parametrize(
    "data, expected_str, expected_tex",
    [
        pytest.param("X(He5)", "He5", "^{5}He", id="isotope_1"),
        pytest.param("X(C12)", "C12", "^{12}C", id="isotope_2"),
        pytest.param("X(A1B2C3)", "AB2C3", "A_{}B_{2}C_{3}", id="molecule_1"),
        pytest.param("X(C12H24)", "C12H24", "C_{12}H_{24}", id="molecule_2"),
        pytest.param("X(H2O)", "H2O", "H_{2}O_{}", id="molecule_2"),
        pytest.param("X(ash)", "ash", "ash", id="descriptive_name"),
    ],
)
def test_Substance_to_str(data, expected_str, expected_tex):
    sub = Substance(data)
    assert str(sub) == expected_str
    assert sub.to_tex() == expected_tex

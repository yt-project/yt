from yt.testing import assert_allclose_units, assert_equal, fake_random_ds
from yt.utilities.chemical_formulas import ChemicalFormula
from yt.utilities.physical_ratios import _primordial_mass_fraction


def test_default_species_fields():

    # Test default case (no species fields)

    ds = fake_random_ds(32)
    sp = ds.sphere("c", (0.2, "unitary"))

    mu = 0.5924489101195808

    assert_allclose_units(mu * sp["index", "ones"], sp["gas", "mean_molecular_weight"])

    assert ("gas", "H_nuclei_density") not in ds.derived_field_list
    assert ("gas", "He_nuclei_density") not in ds.derived_field_list
    assert ("gas", "El_number_density") not in ds.derived_field_list
    assert ("gas", "H_p1_number_density") not in ds.derived_field_list
    assert ("gas", "He_p2_number_density") not in ds.derived_field_list
    assert ("gas", "H_p0_number_density") not in ds.derived_field_list
    assert ("gas", "He_p0_number_density") not in ds.derived_field_list

    # Test fully ionized case
    dsi = fake_random_ds(32, default_species_fields="ionized")
    spi = dsi.sphere("c", (0.2, "unitary"))
    amu_cgs = dsi.units.physical_constants.amu_cgs

    mueinv = 1.0 * _primordial_mass_fraction["H"] / ChemicalFormula("H").weight
    mueinv *= spi["index", "ones"]
    mueinv += 2.0 * _primordial_mass_fraction["He"] / ChemicalFormula("He").weight
    mupinv = _primordial_mass_fraction["H"] / ChemicalFormula("H").weight
    mupinv *= spi["index", "ones"]
    muainv = _primordial_mass_fraction["He"] / ChemicalFormula("He").weight
    muainv *= spi["index", "ones"]
    mueinv2 = spi["gas", "El_number_density"] * amu_cgs / spi["gas", "density"]
    mupinv2 = spi["gas", "H_p1_number_density"] * amu_cgs / spi["gas", "density"]
    muainv2 = spi["gas", "He_p2_number_density"] * amu_cgs / spi["gas", "density"]

    assert_allclose_units(mueinv, mueinv2)
    assert_allclose_units(mupinv, mupinv2)
    assert_allclose_units(muainv, muainv2)

    assert_equal(spi["gas", "H_p1_number_density"], spi["gas", "H_nuclei_density"])
    assert_equal(spi["gas", "He_p2_number_density"], spi["gas", "He_nuclei_density"])

    mu = 0.5924489101195808

    assert_allclose_units(
        mu * spi["index", "ones"], spi["gas", "mean_molecular_weight"]
    )

    # Test fully neutral case

    dsn = fake_random_ds(32, default_species_fields="neutral")
    spn = dsn.sphere("c", (0.2, "unitary"))
    amu_cgs = dsn.units.physical_constants.amu_cgs

    assert ("gas", "El_number_density") not in ds.derived_field_list

    mupinv = _primordial_mass_fraction["H"] / ChemicalFormula("H").weight
    mupinv *= spn["index", "ones"]
    muainv = _primordial_mass_fraction["He"] / ChemicalFormula("He").weight
    muainv *= spn["index", "ones"]
    mupinv2 = spn["gas", "H_p0_number_density"] * amu_cgs / spn["gas", "density"]
    muainv2 = spn["gas", "He_p0_number_density"] * amu_cgs / spn["gas", "density"]

    assert_allclose_units(mueinv, mueinv2)
    assert_allclose_units(mupinv, mupinv2)
    assert_allclose_units(muainv, muainv2)

    assert_equal(spn["gas", "H_p0_number_density"], spn["gas", "H_nuclei_density"])
    assert_equal(spn["gas", "He_p0_number_density"], spn["gas", "He_nuclei_density"])

    mu = 1.2285402715185552

    assert_allclose_units(
        mu * spn["index", "ones"], spn["gas", "mean_molecular_weight"]
    )

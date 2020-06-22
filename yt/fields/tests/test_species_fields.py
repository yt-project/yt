from yt.testing import requires_file, \
    assert_allclose_units, assert_equal
from yt.utilities.answer_testing.framework import \
    data_dir_load
from yt.utilities.physical_ratios import \
    _primordial_mass_fraction
from yt.utilities.chemical_formulas import \
    ChemicalFormula

sloshing = "GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100"

@requires_file(sloshing)
def test_default_species_fields():
    ds = data_dir_load(sloshing)
    sp = ds.sphere("c", (0.2, "unitary"))
    amu_cgs = ds.units.physical_constants.amu_cgs

    mueinv = 1.0*_primordial_mass_fraction["H"] / \
             ChemicalFormula("H").weight
    mueinv *= sp["index","ones"]
    mueinv += 2.0*_primordial_mass_fraction["He"] / \
              ChemicalFormula("He").weight
    mupinv =  _primordial_mass_fraction["H"] / \
              ChemicalFormula("H").weight
    mupinv *= sp["index","ones"]
    muainv =  _primordial_mass_fraction["He"] / \
              ChemicalFormula("He").weight
    muainv *= sp["index","ones"]
    mueinv2 = sp["gas","El_number_density"]*amu_cgs/sp["gas","density"]
    mupinv2 = sp["gas","H_p1_number_density"]*amu_cgs/sp["gas","density"]
    muainv2 = sp["gas","He_p2_number_density"]*amu_cgs/sp["gas","density"]


    assert_allclose_units(mueinv, mueinv2)
    assert_allclose_units(mupinv, mupinv2)
    assert_allclose_units(muainv, muainv2)

    assert_equal(sp["gas","H_p1_number_density"], sp["gas","H_nuclei_density"])
    assert_equal(sp["gas","He_p2_number_density"], sp["gas","He_nuclei_density"])

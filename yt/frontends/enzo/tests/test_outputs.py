"""
Enzo frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.testing import \
    assert_almost_equal, \
    assert_equal, \
    requires_file, \
    units_override_check, \
    assert_array_equal
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load, \
    AnalyticHaloMassFunctionTest, \
    SimulatedHaloMassFunctionTest
from yt.frontends.enzo.api import EnzoDataset

_fields = ("temperature", "density", "velocity_magnitude",
           "velocity_divergence")

two_sphere_test = 'ActiveParticleTwoSphere/DD0011/DD0011'
active_particle_cosmology = 'ActiveParticleCosmology/DD0046/DD0046'
ecp = "enzo_cosmology_plus/DD0046/DD0046"
g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
toro1d = "ToroShockTube/DD0001/data0001"
kh2d = "EnzoKelvinHelmholtz/DD0011/DD0011"

def check_color_conservation(ds):
    species_names = ds.field_info.species_names
    dd = ds.all_data()
    dens_yt = dd["density"].copy()
    # Enumerate our species here
    for s in sorted(species_names):
        if s == "El": continue
        dens_yt -= dd["%s_density" % s]
    dens_yt -= dd["metal_density"]
    delta_yt = np.abs(dens_yt / dd["density"])

    # Now we compare color conservation to Enzo's color conservation
    dd = ds.all_data()
    dens_enzo = dd["Density"].copy()
    for f in sorted(ds.field_list):
        ff = f[1]
        if not ff.endswith("_Density"):
            continue
        start_strings = ["Electron_", "SFR_", "Forming_Stellar_",
                         "Dark_Matter", "Star_Particle_"]
        if any([ff.startswith(ss) for ss in start_strings]):
            continue
        dens_enzo -= dd[f]
    delta_enzo = np.abs(dens_enzo / dd["Density"])
    return assert_almost_equal, delta_yt, delta_enzo

m7 = "DD0010/moving7_0010"
@requires_ds(m7)
def test_moving7():
    ds = data_dir_load(m7)
    yield assert_equal, str(ds), "moving7_0010"
    for test in small_patch_amr(m7, _fields):
        test_moving7.__name__ = test.description
        yield test

@requires_ds(g30, big_data=True)
def test_galaxy0030():
    ds = data_dir_load(g30)
    yield check_color_conservation(ds)
    yield assert_equal, str(ds), "galaxy0030"
    for test in big_patch_amr(ds, _fields):
        test_galaxy0030.__name__ = test.description
        yield test
    assert_equal(ds.particle_type_counts, {'io': 1124453})

@requires_ds(toro1d)
def test_toro1d():
    ds = data_dir_load(toro1d)
    assert_equal(str(ds), 'data0001')
    for test in small_patch_amr(ds, ds.field_list):
        test_toro1d.__name__ = test.description
        yield test

@requires_ds(kh2d)
def test_kh2d():
    ds = data_dir_load(kh2d)
    assert_equal(str(ds), 'DD0011')
    for test in small_patch_amr(ds, ds.field_list):
        test_toro1d.__name__ = test.description
        yield test

@requires_ds(enzotiny)
def test_simulated_halo_mass_function():
    ds = data_dir_load(enzotiny)
    for finder in ["fof", "hop"]:
        yield SimulatedHaloMassFunctionTest(ds, finder)

@requires_ds(enzotiny)
def test_analytic_halo_mass_function():
    ds = data_dir_load(enzotiny)
    for fit in range(1, 6):
        yield AnalyticHaloMassFunctionTest(ds, fit)

@requires_ds(ecp, big_data=True)
def test_ecp():
    ds = data_dir_load(ecp)
    # Now we test our species fields
    yield check_color_conservation(ds)

@requires_file(enzotiny)
def test_units_override():
    for test in units_override_check(enzotiny):
        yield test

@requires_ds(ecp, big_data=True)
def test_nuclei_density_fields():
    ds = data_dir_load(ecp)
    ad = ds.all_data()
    yield assert_array_equal, ad["H_nuclei_density"], \
      (ad["H_number_density"] + ad["H_p1_number_density"])
    yield assert_array_equal, ad["He_nuclei_density"], \
      (ad["He_number_density"] + ad["He_p1_number_density"] +
       ad["He_p2_number_density"])

@requires_file(enzotiny)
def test_EnzoDataset():
    assert isinstance(data_dir_load(enzotiny), EnzoDataset)

@requires_file(two_sphere_test)
@requires_file(active_particle_cosmology)
def test_active_particle_datasets():
    two_sph = data_dir_load(two_sphere_test)
    assert 'AccretingParticle' in two_sph.particle_types_raw
    assert 'io' not in two_sph.particle_types_raw
    assert 'all' in two_sph.particle_types
    assert_equal(len(two_sph.particle_unions), 1)
    pfields = ['GridID', 'creation_time', 'dynamical_time',
               'identifier', 'level', 'metallicity', 'particle_mass']
    pfields += ['particle_position_%s' % d for d in 'xyz']
    pfields += ['particle_velocity_%s' % d for d in 'xyz']

    acc_part_fields = \
        [('AccretingParticle', pf) for pf in ['AccretionRate'] + pfields]

    real_acc_part_fields = sorted(
        [f for f in two_sph.field_list if f[0] == 'AccretingParticle'])
    assert_equal(acc_part_fields, real_acc_part_fields)


    apcos = data_dir_load(active_particle_cosmology)
    assert_equal(['CenOstriker', 'DarkMatter'], apcos.particle_types_raw)
    assert 'all' in apcos.particle_unions

    apcos_fields = [('CenOstriker', pf) for pf in pfields]

    real_apcos_fields = sorted(
        [f for f in apcos.field_list if f[0] == 'CenOstriker'])

    assert_equal(apcos_fields, real_apcos_fields)

    assert_equal(apcos.particle_type_counts,
                 {'CenOstriker': 899755, 'DarkMatter': 32768})

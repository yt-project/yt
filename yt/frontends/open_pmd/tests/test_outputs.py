"""
openPMD frontend tests



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2016, Fabian Koller (HZDR).
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from yt.frontends.open_pmd.data_structures import \
    OpenPMDDataset
from yt.testing import \
    assert_almost_equal, \
    assert_equal, \
    assert_array_equal, \
    requires_file
from yt.utilities.answer_testing.framework import \
    data_dir_load

twoD = "example-2d/hdf5/data00000100.h5"
threeD = "example-3d/hdf5/data00000100.h5"


@requires_file(threeD)
def test_3d_out():
    ds = data_dir_load(threeD)
    field_list = [('all', 'particle_charge'),
                  ('all', 'particle_mass'),
                  ('all', 'particle_momentum_x'),
                  ('all', 'particle_momentum_y'),
                  ('all', 'particle_momentum_z'),
                  ('all', 'particle_positionCoarse_x'),
                  ('all', 'particle_positionCoarse_y'),
                  ('all', 'particle_positionCoarse_z'),
                  ('all', 'particle_positionOffset_x'),
                  ('all', 'particle_positionOffset_y'),
                  ('all', 'particle_positionOffset_z'),
                  ('all', 'particle_weighting'),
                  ('io', 'particle_charge'),
                  ('io', 'particle_mass'),
                  ('io', 'particle_momentum_x'),
                  ('io', 'particle_momentum_y'),
                  ('io', 'particle_momentum_z'),
                  ('io', 'particle_positionCoarse_x'),
                  ('io', 'particle_positionCoarse_y'),
                  ('io', 'particle_positionCoarse_z'),
                  ('io', 'particle_positionOffset_x'),
                  ('io', 'particle_positionOffset_y'),
                  ('io', 'particle_positionOffset_z'),
                  ('io', 'particle_weighting'),
                  ('openPMD', 'E_x'),
                  ('openPMD', 'E_y'),
                  ('openPMD', 'E_z'),
                  ('openPMD', 'rho')]
    domain_dimensions = [26, 26, 201]
    domain_width = [2.08e-05, 2.08e-05, 2.01e-05]

    assert isinstance(ds, OpenPMDDataset)
    yield assert_equal, str(ds), "data00000100.h5"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.particle_types_raw, ('io',)
    assert "all" in ds.particle_unions
    yield assert_array_equal, ds.field_list, field_list
    yield assert_array_equal, ds.domain_domensions, domain_dimensions
    yield assert_almost_equal, ds.current_time, 3.28471214521e-14
    yield assert_almost_equal, ds.domain_right_edge - ds.domain_left_edge, domain_width


@requires_file(twoD)
def test_2d_out():
    ds = data_dir_load(twoD)
    field_list = [('Hydrogen1+', 'particle_charge'),
                  ('Hydrogen1+', 'particle_mass'),
                  ('Hydrogen1+', 'particle_momentum_x'),
                  ('Hydrogen1+', 'particle_momentum_y'),
                  ('Hydrogen1+', 'particle_momentum_z'),
                  ('Hydrogen1+', 'particle_positionCoarse_x'),
                  ('Hydrogen1+', 'particle_positionCoarse_y'),
                  ('Hydrogen1+', 'particle_positionCoarse_z'),
                  ('Hydrogen1+', 'particle_positionOffset_x'),
                  ('Hydrogen1+', 'particle_positionOffset_y'),
                  ('Hydrogen1+', 'particle_positionOffset_z'),
                  ('Hydrogen1+', 'particle_weighting'),
                  ('all', 'particle_charge'),
                  ('all', 'particle_mass'),
                  ('all', 'particle_momentum_x'),
                  ('all', 'particle_momentum_y'),
                  ('all', 'particle_momentum_z'),
                  ('all', 'particle_positionCoarse_x'),
                  ('all', 'particle_positionCoarse_y'),
                  ('all', 'particle_positionCoarse_z'),
                  ('all', 'particle_positionOffset_x'),
                  ('all', 'particle_positionOffset_y'),
                  ('all', 'particle_positionOffset_z'),
                  ('all', 'particle_weighting'),
                  ('electrons', 'particle_charge'),
                  ('electrons', 'particle_mass'),
                  ('electrons', 'particle_momentum_x'),
                  ('electrons', 'particle_momentum_y'),
                  ('electrons', 'particle_momentum_z'),
                  ('electrons', 'particle_positionCoarse_x'),
                  ('electrons', 'particle_positionCoarse_y'),
                  ('electrons', 'particle_positionCoarse_z'),
                  ('electrons', 'particle_positionOffset_x'),
                  ('electrons', 'particle_positionOffset_y'),
                  ('electrons', 'particle_positionOffset_z'),
                  ('electrons', 'particle_weighting'),
                  ('openPMD', 'B_x'),
                  ('openPMD', 'B_y'),
                  ('openPMD', 'B_z'),
                  ('openPMD', 'E_x'),
                  ('openPMD', 'E_y'),
                  ('openPMD', 'E_z'),
                  ('openPMD', 'J_x'),
                  ('openPMD', 'J_y'),
                  ('openPMD', 'J_z'),
                  ('openPMD', 'rho')]
    domain_dimensions = [51, 201, 1]
    domain_width = [3.06e-05, 2.01e-05, 1e+0]

    assert isinstance(ds, OpenPMDDataset)
    yield assert_equal, str(ds), "data00000100.h5"
    yield assert_equal, ds.dimensionality, 2
    yield assert_equal, ds.particle_types_raw, ('Hydrogen1+', 'electrons')
    assert "all" in ds.particle_unions
    yield assert_array_equal, ds.field_list, field_list
    yield assert_array_equal, ds.domain_domensions, domain_dimensions
    yield assert_almost_equal, ds.current_time, 3.29025596712e-14
    yield assert_almost_equal, ds.domain_right_edge - ds.domain_left_edge, domain_width

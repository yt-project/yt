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

from yt.testing import \
    assert_almost_equal, \
    assert_equal, \
    assert_array_equal, \
    requires_file

from yt.utilities.answer_testing.framework import \
    FieldValuesTest, \
    requires_ds, \
    data_dir_load

from yt.frontends.open_pmd.data_structures import \
    OpenPMDDataset

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
    domain_domensions = [26, 26, 201]
    domain_width = [2.08e-05, 2.08e-05, 2.01e-05]

    assert isinstance(ds, OpenPMDDataset)
    yield assert_equal, str(ds), "data00000100.h5"
    yield assert_equal, ds.dimensionality, 3
    yield assert_equal, ds.current_time, 3.28471214521e-14
    yield assert_equal, ds.particle_types_raw, ('io',)
    assert "all" in ds.particle_unions
    yield assert_array_equal, ds.field_list, field_list
    yield assert_array_equal, ds.domain_domensions, domain_domensions
    yield assert_almost_equal, ds.domain_right_edge - ds.domain_left_edge, domain_width
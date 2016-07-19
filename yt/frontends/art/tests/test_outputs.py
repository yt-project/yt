"""
ART frontend tests using D9p a=0.500




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.testing import \
    requires_file, \
    assert_equal, \
    units_override_check, \
    assert_almost_equal
from yt.units.yt_array import \
    YTQuantity
from yt.utilities.answer_testing.framework import \
    requires_ds, \
    FieldValuesTest, \
    PixelizedProjectionValuesTest, \
    data_dir_load
from yt.frontends.art.api import ARTDataset

_fields = (
    ("gas", "density"),
    ("gas", "temperature"),
    ("all", "particle_mass"),
    ("all", "particle_position_x")
)

d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"

@requires_ds(d9p, big_data=True)
def test_d9p():
    ds = data_dir_load(d9p)
    ds.index
    yield assert_equal, str(ds), "10MpcBox_HartGal_csf_a0.500.d"
    dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for dobj_name in dso:
                for weight_field in [None, "density"]:
                    if field[0] not in ds.particle_types:
                        yield PixelizedProjectionValuesTest(
                            d9p, axis, field, weight_field,
                            dobj_name)
            yield FieldValuesTest(d9p, field, dobj_name)

    ad = ds.all_data()
    # 'Ana' variable values output from the ART Fortran 'ANA' analysis code
    AnaNStars = 6255
    yield assert_equal, ad[('stars', 'particle_type')].size, AnaNStars
    yield assert_equal, ad[('specie4', 'particle_type')].size, AnaNStars

    # The *real* asnwer is 2833405, but yt misses one particle since it lives
    # on a domain boundary. See issue 814. When that is fixed, this test
    # will need to be updated
    AnaNDM = 2833404
    yield assert_equal, ad[('darkmatter', 'particle_type')].size, AnaNDM
    yield assert_equal, (ad[('specie0', 'particle_type')].size +
                         ad[('specie1', 'particle_type')].size +
                         ad[('specie2', 'particle_type')].size +
                         ad[('specie3', 'particle_type')].size), AnaNDM

    for spnum in range(5):
        npart_read = ad['specie%s' % spnum, 'particle_type'].size
        npart_header = ds.particle_type_counts['specie%s' % spnum]
        if spnum == 3:
            # see issue 814
            npart_read += 1
        assert_equal(npart_read, npart_header)

    AnaBoxSize = YTQuantity(7.1442196564, 'Mpc')
    AnaVolume = YTQuantity(364.640074656, 'Mpc**3')
    Volume = 1
    for i in ds.domain_width.in_units('Mpc'):
        yield assert_almost_equal, i, AnaBoxSize
        Volume *= i
    yield assert_almost_equal, Volume, AnaVolume

    AnaNCells = 4087490
    yield assert_equal, len(ad[('index', 'cell_volume')]), AnaNCells

    AnaTotDMMass = YTQuantity(1.01191786808255e+14, 'Msun')
    yield (assert_almost_equal,
           ad[('darkmatter', 'particle_mass')].sum().in_units('Msun'),
           AnaTotDMMass)

    AnaTotStarMass = YTQuantity(1776701.3990607238, 'Msun')
    yield (assert_almost_equal,
           ad[('stars', 'particle_mass')].sum().in_units('Msun'),
           AnaTotStarMass)

    AnaTotStarMassInitial = YTQuantity(2423468.2801332865, 'Msun')
    yield (assert_almost_equal,
           ad[('stars', 'particle_mass_initial')].sum().in_units('Msun'),
           AnaTotStarMassInitial)

    AnaTotGasMass = YTQuantity(1.7826982029216785e+13, 'Msun')
    yield (assert_almost_equal, ad[('gas', 'cell_mass')].sum().in_units('Msun'),
           AnaTotGasMass)

    AnaTotTemp = YTQuantity(150219844793.39072, 'K')  # just leaves
    yield assert_equal, ad[('gas', 'temperature')].sum(), AnaTotTemp


@requires_file(d9p)
def test_ARTDataset():
    assert isinstance(data_dir_load(d9p), ARTDataset)

@requires_file(d9p)
def test_units_override():
    for test in units_override_check(d9p):
        yield test

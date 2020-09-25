from yt.frontends.art.api import ARTDataset
from yt.testing import (
    ParticleSelectionComparison,
    assert_almost_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.units.yt_array import YTQuantity
from yt.utilities.answer_testing.framework import (
    FieldValuesTest,
    PixelizedProjectionValuesTest,
    data_dir_load,
    requires_ds,
)

_fields = (
    ("gas", "density"),
    ("gas", "temperature"),
    ("all", "particle_mass"),
    ("all", "particle_position_x"),
)

d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"
dmonly = "DMonly/PMcrs0.0100.DAT"


@requires_ds(d9p, big_data=True)
def test_d9p():
    ds = data_dir_load(d9p)
    ds.index
    assert_equal(str(ds), "10MpcBox_HartGal_csf_a0.500.d")
    dso = [None, ("sphere", ("max", (0.1, "unitary")))]
    for field in _fields:
        for axis in [0, 1, 2]:
            for dobj_name in dso:
                for weight_field in [None, ("gas", "density")]:
                    if field[0] not in ds.particle_types:
                        yield PixelizedProjectionValuesTest(
                            d9p, axis, field, weight_field, dobj_name
                        )
            yield FieldValuesTest(d9p, field, dobj_name)


@requires_file(d9p)
def test_d9p_ana_values():
    ds = data_dir_load(d9p)
    ad = ds.all_data()
    # 'Ana' variable values output from the ART Fortran 'ANA' analysis code
    AnaNStars = 6255
    assert_equal(ad[("stars", "particle_type")].size, AnaNStars)
    assert_equal(ad[("specie4", "particle_type")].size, AnaNStars)

    # The *real* asnwer is 2833405, but yt misses one particle since it lives
    # on a domain boundary. See issue 814. When that is fixed, this test
    # will need to be updated
    AnaNDM = 2_833_404
    assert_equal(ad[("darkmatter", "particle_type")].size, AnaNDM)
    assert_equal(
        (
            ad[("specie0", "particle_type")].size
            + ad[("specie1", "particle_type")].size
            + ad[("specie2", "particle_type")].size
            + ad[("specie3", "particle_type")].size
        ),
        AnaNDM,
    )

    for spnum in range(5):
        npart_read = ad[f"specie{spnum}", "particle_type"].size
        npart_header = ds.particle_type_counts[f"specie{spnum}"]
        if spnum == 3:
            # see issue 814
            npart_read += 1
        assert_equal(npart_read, npart_header)

    AnaBoxSize = YTQuantity(7.144_219_656_4, "Mpc")
    AnaVolume = YTQuantity(364.640_074_656, "Mpc**3")
    Volume = 1
    for i in ds.domain_width.in_units("Mpc"):
        assert_almost_equal(i, AnaBoxSize)
        Volume *= i
    assert_almost_equal(Volume, AnaVolume)

    AnaNCells = 4_087_490
    assert_equal(len(ad[("index", "cell_volume")]), AnaNCells)

    AnaTotDMMass = YTQuantity(1.011_917_868_082_55e14, "Msun")
    assert_almost_equal(
        ad[("darkmatter", "particle_mass")].sum().in_units("Msun"), AnaTotDMMass
    )

    AnaTotStarMass = YTQuantity(1_776_701.399_060_723_8, "Msun")
    assert_almost_equal(
        ad[("stars", "particle_mass")].sum().in_units("Msun"), AnaTotStarMass
    )

    AnaTotStarMassInitial = YTQuantity(2_423_468.280_133_286_5, "Msun")
    assert_almost_equal(
        ad[("stars", "particle_mass_initial")].sum().in_units("Msun"),
        AnaTotStarMassInitial,
    )

    AnaTotGasMass = YTQuantity(1.782_698_202_921_678_5e13, "Msun")
    assert_almost_equal(ad[("gas", "cell_mass")].sum().in_units("Msun"), AnaTotGasMass)

    AnaTotTemp = YTQuantity(150_219_844_793.39072, "K")  # just leaves
    assert_almost_equal(ad[("gas", "temperature")].sum(), AnaTotTemp, decimal=4)


@requires_file(d9p)
def test_ARTDataset():
    assert isinstance(data_dir_load(d9p), ARTDataset)


@requires_file(d9p)
def test_units_override():
    units_override_check(d9p)


@requires_file(dmonly)
def test_particle_selection():
    ds = data_dir_load(dmonly)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()

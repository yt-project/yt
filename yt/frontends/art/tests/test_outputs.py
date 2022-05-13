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
    ("all", "particle_velocity_y"),
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
        for axis in [0, 1]:
            for dobj_name in dso:
                for weight_field in [None, ("gas", "density")]:
                    if field[0] not in ds.particle_types:
                        yield PixelizedProjectionValuesTest(
                            d9p, axis, field, weight_field, dobj_name
                        )
                yield FieldValuesTest(
                    d9p,
                    field,
                    obj_type=dobj_name,
                    particle_type=field[0] in ds.particle_types,
                )


@requires_ds(d9p, big_data=True)
def test_d9p_global_values():
    ds = data_dir_load(d9p)
    ad = ds.all_data()
    # 'Ana' variable values output from the ART Fortran 'ANA' analysis code
    AnaNStars = 6255
    assert_equal(ad[("stars", "particle_type")].size, AnaNStars)
    assert_equal(ad[("specie4", "particle_type")].size, AnaNStars)

    # The *real* answer is 2833405, but yt misses one particle since it lives
    # on a domain boundary. See issue 814. When that is fixed, this test
    # will need to be updated
    AnaNDM = 2833404
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

    AnaBoxSize = YTQuantity(7.1442196564, "Mpc")
    AnaVolume = YTQuantity(364.640074656, "Mpc**3")
    Volume = 1
    for i in ds.domain_width.in_units("Mpc"):
        assert_almost_equal(i, AnaBoxSize)
        Volume *= i
    assert_almost_equal(Volume, AnaVolume)

    AnaNCells = 4087490
    assert_equal(len(ad[("index", "cell_volume")]), AnaNCells)

    AnaTotDMMass = YTQuantity(1.01191786808255e14, "Msun")
    assert_almost_equal(
        ad[("darkmatter", "particle_mass")].sum().in_units("Msun"), AnaTotDMMass
    )

    AnaTotStarMass = YTQuantity(1776701.3990607238, "Msun")
    assert_almost_equal(
        ad[("stars", "particle_mass")].sum().in_units("Msun"), AnaTotStarMass
    )

    AnaTotStarMassInitial = YTQuantity(2423468.2801332865, "Msun")
    assert_almost_equal(
        ad[("stars", "particle_mass_initial")].sum().in_units("Msun"),
        AnaTotStarMassInitial,
    )

    AnaTotGasMass = YTQuantity(1.7826982029216785e13, "Msun")
    assert_almost_equal(ad[("gas", "cell_mass")].sum().in_units("Msun"), AnaTotGasMass)

    AnaTotTemp = YTQuantity(150219844793.3907, "K")  # just leaves
    assert_almost_equal(ad[("gas", "temperature")].sum().in_units("K"), AnaTotTemp)


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

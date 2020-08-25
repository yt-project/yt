import numpy as np
import pytest

from yt.frontends.art.api import ARTDataset
from yt.testing import (
    ParticleSelectionComparison,
    assert_almost_equal,
    assert_equal,
    requires_file,
    units_override_check,
)
from yt.units.yt_array import YTQuantity
from yt.utilities.answer_testing.answer_tests import (
    field_values,
    pixelized_projection_values,
)

# Test data
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"

a_list = [0, 1, 2]
d_list = [None, ("sphere", ("max", (0.1, "unitary")))]
w_list = [None, "density"]
f_list = [
    ("gas", "density"),
    ("gas", "temperature"),
    ("all", "particle_mass"),
    ("all", "particle_position_x"),
]


@pytest.mark.answer_test
class TestArt:
    answer_file = None
    saved_hashes = None

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @pytest.mark.parametrize("ds", [d9p], indirect=True)
    @pytest.mark.parametrize("f", f_list, indirect=True)
    @pytest.mark.parametrize("w", w_list, indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    @pytest.mark.parametrize("a", a_list, indirect=True)
    def test_d9p_ppv(self, f, d, a, w, ds):
        ds.index
        particle_type = f[0] in ds.particle_types
        if not particle_type:
            ppv = pixelized_projection_values(ds, a, f, w, d)
            self.hashes.update({"pixelized_projection_values": ppv})
        # So we have something to save for this test in the answer file
        else:
            self.hashes.update({"pixelized_projection_values": np.array(-1)})

    @pytest.mark.parametrize("ds", [d9p], indirect=True)
    @pytest.mark.parametrize("f", f_list, indirect=True)
    @pytest.mark.parametrize("d", d_list, indirect=True)
    def test_d9p_fv(self, f, d, ds):
        ds.index
        # Does it need the particle_type argument passed?
        fv = field_values(ds, f, d)
        self.hashes.update({"field_values": fv})

    @pytest.mark.big_data
    @pytest.mark.parametrize("ds", [d9p], indirect=True)
    def test_d9p_no_params(self, ds):
        ds.index
        ad = ds.all_data()
        # 'Ana' variable values output from the ART Fortran 'ANA' analysis code
        AnaNStars = 6255
        assert_equal(ad[("stars", "particle_type")].size, AnaNStars)
        assert_equal(ad[("specie4", "particle_type")].size, AnaNStars)
        # The *real* asnwer is 2833405, but yt misses one particle since it lives
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
        assert_almost_equal(
            ad[("gas", "cell_mass")].sum().in_units("Msun"), AnaTotGasMass
        )
        AnaTotTemp = YTQuantity(150219844793.39072, "K")  # just leaves
        assert_equal(ad[("gas", "temperature")].sum(), AnaTotTemp)

    @pytest.mark.parametrize("ds", [d9p], indirect=True)
    def test_ARTDataset(self, ds):
        assert isinstance(ds, ARTDataset)

    @requires_file(d9p)
    def test_units_override(self):
        units_override_check(d9p)

    @pytest.mark.parametrize("ds", [d9p], indirect=True)
    def test_particle_selection(self, ds):
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

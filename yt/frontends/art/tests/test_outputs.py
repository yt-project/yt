"""
Title: test_art.py
Purpose: ART frontend tests using D9p a=0.500
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this software.
"""
from yt.testing import \
    units_override_check, \
    assert_almost_equal
from yt.units.yt_array import \
    YTQuantity
from yt.frontends.art.api import ARTDataset

import framework as fw
import utils


# Test data
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


#============================================
#                   TestArt
#============================================
class TestArt(fw.AnswerTest):
    """
    Container for art frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_d9p
    #-----
    def test_d9p(self):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        fields = (
            ("gas", "density"),
            ("gas", "temperature"),
            ("all", "particle_mass"),
            ("all", "particle_position_x")
        )
        ppv_hd = b''
        fv_hd = b''
        ds = utils.data_dir_load(d9p)
        ds.index
        assert str(ds) == "10MpcBox_HartGal_csf_a0.500.d"
        dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
        for field in fields:
            for axis in [0, 1, 2]:
                for dobj_name in dso:
                    for weight_field in [None, "density"]:
                        if field[0] not in ds.particle_types:
                            ppv_hd += self.pixelized_projection_values_test(
                                ds, axis, field, weight_field,
                                dobj_name)
                if field[0] == "all":
                    particle_type = True
                else:
                    particle_type = False
                fv_hd += self.field_values_test(ds, field, dobj_name, particle_type=particle_type)
        ad = ds.all_data()
        # 'Ana' variable values output from the ART Fortran 'ANA' analysis code
        AnaNStars = 6255
        assert ad[('stars', 'particle_type')].size == AnaNStars
        assert ad[('specie4', 'particle_type')].size == AnaNStars
        # The *real* asnwer is 2833405, but yt misses one particle since it lives
        # on a domain boundary. See issue 814. When that is fixed, this test
        # will need to be updated
        AnaNDM = 2833404
        assert ad[('darkmatter', 'particle_type')].size == AnaNDM
        assert (ad[('specie0', 'particle_type')].size +\
                  ad[('specie1', 'particle_type')].size +\
                  ad[('specie2', 'particle_type')].size +\
                  ad[('specie3', 'particle_type')].size) == AnaNDM
        for spnum in range(5):
            npart_read = ad['specie%s' % spnum, 'particle_type'].size
            npart_header = ds.particle_type_counts['specie%s' % spnum]
            if spnum == 3:
                # see issue 814
                npart_read += 1
            assert npart_read == npart_header
        AnaBoxSize = YTQuantity(7.1442196564, 'Mpc')
        AnaVolume = YTQuantity(364.640074656, 'Mpc**3')
        Volume = 1
        for i in ds.domain_width.in_units('Mpc'):
            assert_almost_equal(i, AnaBoxSize)
            Volume *= i
        assert_almost_equal(Volume, AnaVolume)
        AnaNCells = 4087490
        assert len(ad[('index', 'cell_volume')]) == AnaNCells
        AnaTotDMMass = YTQuantity(1.01191786808255e+14, 'Msun')
        assert_almost_equal(
            ad[('darkmatter', 'particle_mass')].sum().in_units('Msun'),
            AnaTotDMMass)
        AnaTotStarMass = YTQuantity(1776701.3990607238, 'Msun')
        assert_almost_equal(ad[('stars', 'particle_mass')].sum().in_units('Msun'),
                            AnaTotStarMass)
        AnaTotStarMassInitial = YTQuantity(2423468.2801332865, 'Msun')
        assert_almost_equal(
            ad[('stars', 'particle_mass_initial')].sum().in_units('Msun'),
            AnaTotStarMassInitial)
        AnaTotGasMass = YTQuantity(1.7826982029216785e+13, 'Msun')
        assert_almost_equal(ad[('gas', 'cell_mass')].sum().in_units('Msun'),
                            AnaTotGasMass)
        AnaTotTemp = YTQuantity(150219844793.39072, 'K')  # just leaves
        assert ad[('gas', 'temperature')].sum() == AnaTotTemp
        hashes = {'pixelized_projection_values' : utils.generate_hash(ppv_hd),
            'field_values' : utils.generate_hash(fv_hd)
        }
        utils.handle_hashes(self.save_dir, 'art', hashes, self.answer_store)

    #-----
    # test_ARTDataset
    #-----
    def test_ARTDataset(self):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        assert isinstance(utils.data_dir_load(d9p), ARTDataset)

    #-----
    # test_units_override
    #-----
    def test_units_override(self):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        units_override_check(d9p)

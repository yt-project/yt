"""
Title: test_artio.py
Purpose: Contains ARTIO frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from yt.convenience import load
from yt.frontends.artio.api import ARTIODataset
from yt.testing import \
    assert_allclose_units, \
    assert_equal, \
    requires_file, \
    units_override_check
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Data file
sizmbhloz = "sizmbhloz-clref04SNth-rs9_a0.9011/"
sizmbhloz += "sizmbhloz-clref04SNth-rs9_a0.9011.art"


#============================================
#                 TestArtIo
#============================================
class TestArtIo(fw.AnswerTest):
    """
    Container for ARTIO answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_sizmbhloz
    #-----
    @utils.requires_ds(sizmbhloz)
    def test_sizmbhloz(self, ds_sizmbhloz):
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
        # Load data
        ds_sizmbhloz.max_range = 1024*1024
        # Set up test parameters
        dso = [None, ("sphere", ("max", (0.1, 'unitary')))]
        axes = [0, 1, 2]
        weight_fields = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
                   ("deposit", "all_density"), ("deposit", "all_count"))
        # Set up hex digests
        ppv_hd = b''
        fv_hd = b''
        # Run tests
        for dobj_name in dso:
            for field in fields:
                for axis in axes:
                    for weight_field in weight_fields:
                        ppv_hd += self.pixelized_projection_values_test(
                            ds_sizmbhloz, axis, field, weight_field,
                            dobj_name
                        )
                fv_hd += self.field_values_test(ds_sizmbhloz, field, dobj_name)
            dobj = utils.create_obj(ds_sizmbhloz, dobj_name)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert_equal(s1, s2)
        assert_equal(ds_sizmbhloz.particle_type_counts, {'N-BODY': 100000, 'STAR': 110650})
        # Save or compare hashes
        hashes = {}
        hashes['pixelized_projection_values'] = utils.generate_hash(ppv_hd)
        hashes['field_values'] = utils.generate_hash(fv_hd)
        utils.handle_hashes(self.save_dir, 'artio', hashes, self.answer_store)

    #-----
    # test_ARTIODataset
    #-----
    @requires_file(sizmbhloz)
    def test_ARTIODataset(self, ds_sizmbhloz):
        """
        Makes sure the loaded data is the proper type.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns
        -------
            pass
        """
        assert isinstance(ds_sizmbhloz, ARTIODataset)

    #-----
    # test_units_override
    #-----
    @requires_file(sizmbhloz)
    def test_units_override(self, ds_sizmbhloz):
        """
        Performs the units override test.

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
        units_override_check(ds_sizmbhloz, sizmbhloz)

    #-----
    # test_particle_derived_field
    #-----
    @requires_file(sizmbhloz)
    def test_particle_derived_field(self, ds_sizmbhloz):
        """
        Defines a derived field and makes sure that it works.

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
        def star_age_alias(field, data):
            # test to make sure we get back data in the correct units
            # during field detection
            return data['STAR', 'age'].in_units('Myr')

        ds_sizmbhloz.add_field(("STAR", "new_field"), function=star_age_alias,
                     units='Myr', sampling_type="particle")
        ad = ds_sizmbhloz.all_data()
        assert_allclose_units(ad['STAR', 'age'].in_units("Myr"),
                              ad["STAR", "new_field"])

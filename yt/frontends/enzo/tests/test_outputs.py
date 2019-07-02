"""
Title:   test_enzo.py
Purpose: Contains Enzo frontend tests
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import yt
from yt.frontends.enzo.api import EnzoDataset
from yt.utilities.answer_testing import framework as fw
from yt.utilities.answer_testing import utils


# Files containing data to be used in tests. Paths are relative to
# yt test_data_dir
toro1d = "ToroShockTube/DD0001/data0001"
kh2d = "EnzoKelvinHelmholtz/DD0011/DD0011"
m7 = "DD0010/moving7_0010"
g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"
enzotiny = "enzo_tiny_cosmology/DD0046/DD0046"
ecp = "enzo_cosmology_plus/DD0046/DD0046"


#============================================
#                 TestEnzo
#============================================
class TestEnzo(fw.AnswerTest):
    """
    Container for Enzo answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_toro1d
    #-----
    #@utils.requires_ds(toro1d)
    def test_toro1d(self):
        # Load data
        ds = utils.data_dir_load(toro1d)
        # Make sure we have the right data file
        assert str(ds) == 'data0001'
        # Set up arrays for testing
        axes = [0, 1, 2]
        input_center = "max"
        ds_objs = [None, ("sphere", (input_center, (0.1, 'unitary')))]
        weight_fields = [None, "density"]
        fields = ds.field_list
        # Grid hierarchy test
        gh_hd = utils.generate_hash(
            self.grid_hierarchy_test(ds)
        )
        # Parentage relationships test
        pr_hd = utils.generate_hash(
            self.parentage_relationships_test(ds)
        )
        # Grid values, projection values, and field values tests
        # For these tests it might be possible for a frontend to not
        # run them for every field, axis, or data source, so the tests
        # are written to work with one instance of each of those. These
        # loops therefore combine the results of each desired combo.
        # This mirrors the way that it was done in the original test
        gv_hd = b''
        pv_hd = b''
        fv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weight_fields:
                        pv_hd += self.projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
                    fv_hd += self.field_values_test(ds,
                        field,
                        dobj_name
                    )
        # Hash the final byte arrays
        gv_hd = utils.generate_hash(gv_hd)
        pv_hd = utils.generate_hash(pv_hd)
        fv_hd = utils.generate_hash(fv_hd)
        # Prepare hashes for either writing or comparing
        hashes = {}
        hashes['grid_hierarchy'] = gh_hd
        hashes['parentage_relationships'] = pr_hd
        hashes['grid_values'] = gv_hd
        hashes['projection_values'] = pv_hd
        hashes['field_values'] = fv_hd
        # Save answer
        if self.answer_store:
            utils.store_hashes('toro1d', hashes)
        # Compare to already saved answer
        else:
            saved_hashes = utils.load_hashes('toro1d')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_kh2d
    #-----
    #@requires_ds(kh2d)
    def test_kh2d(self):
        # Load data
        ds = utils.data_dir_load(kh2d)
        # Make sure we're dealing with the right file
        assert str(ds) == 'DD0011'
        # Set up arrays for testing
        axes = [0, 1, 2]
        input_center = "max"
        ds_objs = [None, ("sphere", (input_center, (0.1, 'unitary')))]
        weight_fields = [None, "density"]
        fields = ds.field_list
        # Grid hierarchy test
        gh_hd = utils.generate_hash(
            self.grid_hierarchy_test(ds)
        )
        # Parentage relationships test
        pr_hd = utils.generate_hash(
            self.parentage_relationships_test(ds)
        )
        # Grid values, projection values, and field values tests
        # For these tests it might be possible for a frontend to not
        # run them for every field, axis, or data source, so the tests
        # are written to work with one instance of each of those. These
        # loops therefore combine the results of each desired combo.
        # This mirrors the way that it was done in the original test
        gv_hd = b''
        pv_hd = b''
        fv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weight_fields:
                        pv_hd += self.projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
                    fv_hd += self.field_values_test(ds,
                        field,
                        dobj_name
                    )
        # Hash the final byte arrays
        gv_hd = utils.generate_hash(gv_hd)
        pv_hd = utils.generate_hash(pv_hd)
        fv_hd = utils.generate_hash(fv_hd)
        # Prepare hashes for either writing or comparing
        hashes = {}
        hashes['grid_hierarchy'] = gh_hd
        hashes['parentage_relationships'] = pr_hd
        hashes['grid_values'] = gv_hd
        hashes['projection_values'] = pv_hd
        hashes['field_values'] = fv_hd
        # Save answer
        if self.answer_store:
            utils.store_hashes('kh2d', hashes)
        # Compare to already saved answer
        else:
            saved_hashes = utils.load_hashes('kh2d')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_moving7
    #-----
    #@utils.requires_ds(m7)
    def test_moving7(self):
        # Load data
        ds = utils.data_dir_load(m7)
        # Make sure we have the right data file
        assert str(ds) == 'moving7_0010'
        # Set up arrays for testing
        axes = [0, 1, 2]
        input_center = "max"
        ds_objs = [None, ("sphere", (input_center, (0.1, 'unitary')))]
        weight_fields = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
            "velocity_divergence"
        )
        # Grid hierarchy test
        gh_hd = utils.generate_hash(
            self.grid_hierarchy_test(ds)
        )
        # Parentage relationships test
        pr_hd = utils.generate_hash(
            self.parentage_relationships_test(ds)
        )
        # Grid values, projection values, and field values tests
        # For these tests it might be possible for a frontend to not
        # run them for every field, axis, or data source, so the tests
        # are written to work with one instance of each of those. These
        # loops therefore combine the results of each desired combo.
        # This mirrors the way that it was done in the original test
        gv_hd = b''
        pv_hd = b''
        fv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weight_fields:
                        pv_hd += self.projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
                    fv_hd += self.field_values_test(ds,
                        field,
                        dobj_name
                    )
        # Hash the final byte arrays
        gv_hd = utils.generate_hash(gv_hd)
        pv_hd = utils.generate_hash(pv_hd)
        fv_hd = utils.generate_hash(fv_hd)
        # Prepare hashes for either writing or comparing
        hashes = {}
        hashes['grid_hierarchy'] = gh_hd
        hashes['parentage_relationships'] = pr_hd
        hashes['grid_values'] = gv_hd
        #hashes['projection_values'] = pv_hd
        #hashes['field_values'] = fv_hd
        # Save answer
        if self.answer_store:
            utils.store_hashes('moving7', hashes)
        # Compare to already saved answer
        else:
            saved_hashes = utils.load_hashes('moving7')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_galaxy0030
    #-----
    #@requres_ds(g30, big_data=True)
    def test_galaxy0030(self):
        """
        Tests the galaxy_0030 big data simulation.

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
        ds = utils.data_dir_load(g30)
        # Make sure we have the right data file
        assert str(ds) == 'galaxy0030'
        # Set up arrays for testing
        axes = [0, 1, 2]
        input_center = "max"
        ds_objs = [None, ("sphere", (input_center, (0.1, 'unitary')))]
        weight_fields = [None, "density"]
        fields = ("temperature", "density", "velocity_magnitude",
                   "velocity_divergence")
        # Color conservation test
        self.color_conservation_test(ds)
        # Grid hierarchy test
        gh_hd = utils.generate_hash(
            self.grid_hierarchy_test(ds)
        )
        # Parentage relationships test
        pr_hd = utils.generate_hash(
            self.parentage_relationships_test(ds)
        )
        # Grid values and pixelized projection values tests
        # For these tests it might be possible for a frontend to not
        # run them for every field, axis, or data source, so the tests
        # are written to work with one instance of each of those. These
        # loops therefore combine the results of each desired combo.
        # This mirrors the way that it was done in the original test
        gv_hd = b''
        ppv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weight_fields:
                        ppv_hd += self.pixelized_projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
        # Hash the final byte arrays
        gv_hd = utils.generate_hash(gv_hd)
        ppv_hd = utils.generate_hash(ppv_hd)
        # Prepare hashes for either writing or comparing
        hashes = {}
        hashes['grid_hierarchy'] = gh_hd
        hashes['parentage_relationships'] = pr_hd
        hashes['grid_values'] = gv_hd
        hashes['pixelized_projection_values'] = ppv_hd
        # Save answer
        if self.answer_store:
            utils.store_hashes('galaxy0030', hashes)
        # Compare to already saved answer
        else:
            saved_hashes = utils.load_hashes('galaxy0030')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_simulated_halo_mass_function
    #-----
    #@requires_ds(enzotiny)
    def test_simulated_halo_mass_function(self):
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
        # Load in the data
        ds = utils.data_dir_load(enzotiny)
        # Set up hex digest
        shmf_hd = b''
        # Loop over the different finders
        for finder in ["fof", "hop"]:
            shmf_hd += self.simulated_halo_mass_function_test(ds, finder)
        hashes = {}
        hashes['simulated_halo_mass_function'] = utils.generate_hash(shmf_hd)
        if self.answer_store:
            utils.store_hashes('shmf-enzotiny', hashes)
        else:
            saved_hashes = utils.load_hashes('shmf-enzotiny')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_analytic_halo_mass_function
    #-----
    #@requires_ds(enzotiny)
    def test_analytic_halo_mass_function(self):
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
        # Load in the data
        ds = utils.data_dir_load(enzotiny)
        # Set up hex digest
        ahmf_hd = b''
        # Loop over the different finders
        for fit in range(1,6):
            ahmf_hd += self.analytic_halo_mass_function_test(ds, fit)
        hashes = {}
        hashes['analytic_halo_mass_function'] = utils.generate_hash(ahmf_hd)
        if self.answer_store:
            utils.store_hashes('ahmf-enzotiny', hashes)
        else:
            saved_hashes = utils.load_hashes('ahmf-enzotiny')
            for key, value in hashes.items():
                assert saved_hashes[key] == hashes[key]

    #-----
    # test_ecp
    #-----
    #@requires_ds(ecp, big_data=True)
    def test_ecp(self):
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
        ds = utils.data_dir_load(ecp)
        self.color_conservation_test(ds)

    #-----
    # test_units_override
    #-----
    #@requires_file(enzotiny)
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
        yt.testing.units_override_check(enzotiny)

    #-----
    # test_nuclei_density_fields
    #-----
    #@requires_ds(ecp, big_data=True)
    def test_nuclei_density_fields(self):
        ds = utils.data_dir_load(ecp)
        ad = ds.all_data()
        # Hydrogen
        hd1 = utils.generate_hash(ad["H_nuclei_density"].tostring())
        hd2 = utils.generate_hash((ad["H_number_density"] +
            ad["H_p1_number_density"]).tostring())
        assert hd1 == hd2
        hd1 = utils.generate_hash(ad["He_nuclei_density"].tostring())
        hd2 = utils.generate_hash((ad["He_number_density"] +
            ad["He_p1_number_density"] +
            ad["He_p2_number_density"]).tostring()
        )
        assert hd1 == hd2

    #-----
    # test_enzo_dataset
    #-----
    #@requires_file(enzotiny)
    def test_enzo_dataset(self):
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
        assert isinstance(utils.data_dir_load(enzotiny), EnzoDataset)

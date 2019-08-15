"""
Title: test_gadget_fof.py
Purpose: GadgetFOF frontend tests using gadget_fof datasets
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np

from yt.frontends.gadget_fof.api import GadgetFOFDataset

import framework as fw
import utils

# Test data
g5 = "gadget_fof_halos/groups_005/fof_subhalo_tab_005.0.hdf5"
g42 = "gadget_fof_halos/groups_042/fof_subhalo_tab_042.0.hdf5"
g298 = "gadget_halos/data/groups_298/fof_subhalo_tab_298.0.hdf5"
g56 = "gadget_halos/data/groups_056/fof_subhalo_tab_056.0.hdf5"
g76 = "gadget_halos/data/groups_076/fof_subhalo_tab_076.0.hdf5"


#============================================
#               TestGadgetFOF
#============================================
class TestGadgetFOF(fw.AnswerTest):
    """
    Container for gadget fof frontend tests.

    Attritubes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_fields_g5
    #-----
    def test_fields_g5(self):
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
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_velocity_x",
                   "particle_velocity_y", "particle_velocity_z",
                   "particle_mass", "particle_identifier")
        fv_hd = b''
        ds = utils.data_dir_load(g5)
        for field in fields:
            fv_hd += self.field_values_test(ds, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'g5-fields', hashes, self.answer_store)

    #-----
    # test_fields_g42
    #-----
    def test_fields_g42(self):
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
        fields = ("particle_position_x", "particle_position_y",
                   "particle_position_z", "particle_velocity_x",
                   "particle_velocity_y", "particle_velocity_z",
                   "particle_mass", "particle_identifier")
        fv_hd = b''
        ds = utils.data_dir_load(g42)
        for field in fields:
            fv_hd += self.field_values_test(ds, field, particle_type=True)
        hashes = {'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'g42-fields', hashes, self.answer_store)

    #-----
    # test_GadgetFOFDataset
    #-----
    def test_GadgetFOFDataset(self):
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
        assert isinstance(utils.data_dir_load(g42), GadgetFOFDataset)

    #-----
    # test_subhalos
    #-----
    def test_subhalos(self):
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
        ds = utils.data_dir_load(g298)
        total_sub = 0
        total_int = 0
        for hid in range(0, ds.index.particle_count["Group"]):
            my_h = ds.halo("Group", hid)
            h_ids = my_h["ID"]
            for sid in range(int(my_h["subhalo_number"][0])):
                my_s = ds.halo("Subhalo", (my_h.particle_identifier, sid))
                total_sub += my_s["ID"].size
                total_int += np.intersect1d(h_ids, my_s["ID"]).size
        assert total_sub == total_int

    #-----
    # test_halo_masses
    #-----
    def test_halo_masses(self):
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
        ds = utils.data_dir_load(g298)
        ad = ds.all_data()
        for ptype in ["Group", "Subhalo"]:
            nhalos = ds.index.particle_count[ptype]
            mass = ds.arr(np.zeros(nhalos), "code_mass")
            for i in range(nhalos):
                halo = ds.halo(ptype, i)
                mass[i] = halo.mass
            np.testing.assert_array_equal(ad[ptype, "particle_mass"], mass)

    #-----
    # test_unbalanced_dataset
    #-----
    def test_unbalanced_dataset(self):
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
        ds = utils.data_dir_load(g56)
        halo = ds.halo("Group", 0)
        halo["member_ids"]
        assert True

    #-----
    # test_3file_halo
    #-----
    def test_3file_halo(self):
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
        ds = utils.data_dir_load(g76)
        # this halo's particles are distributed over 3 files with the
        # middle file being empty
        halo = ds.halo("Group", 6)
        halo["member_ids"]
        assert True

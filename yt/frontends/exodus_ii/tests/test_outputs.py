"""
Title: test_exodusii.py
Purpose: Exodus II frontend tests
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import framework as fw
import utils

# Test data
out = "ExodusII/out.e"
out_s002 = "ExodusII/out.e-s002"
gold = "ExodusII/gold.e"
big_data = "MOOSE_sample_data/mps_out.e"


#============================================
#               TestExodusII
#============================================
class TestExodusII(fw.AnswerTest):
    """
    Container for exodus ii frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_out
    #-----
    def test_out(self):
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
        ds = utils.data_dir_load(out)
        field_list = [('all', 'conv_indicator'),
                      ('all', 'conv_marker'),
                      ('all', 'convected'),
                      ('all', 'diffused'),
                      ('connect1', 'conv_indicator'),
                      ('connect1', 'conv_marker'),
                      ('connect1', 'convected'),
                      ('connect1', 'diffused'),
                      ('connect2', 'conv_indicator'),
                      ('connect2', 'conv_marker'),
                      ('connect2', 'convected'),
                      ('connect2', 'diffused')]
        assert str(ds) == "out.e"
        assert ds.dimensionality == 3
        assert ds.current_time == 0.0
        assert ds.parameters['nod_names'] == ['convected', 'diffused']
        assert ds.parameters['num_meshes'] == 2
        assert ds.field_list == field_list

    #-----
    # test_out002
    #-----
    def test_out002(self):
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
        ds = utils.data_dir_load(out_s002)
        field_list = [('all', 'conv_indicator'),
                      ('all', 'conv_marker'),
                      ('all', 'convected'),
                      ('all', 'diffused'),
                      ('connect1', 'conv_indicator'),
                      ('connect1', 'conv_marker'),
                      ('connect1', 'convected'),
                      ('connect1', 'diffused'),
                      ('connect2', 'conv_indicator'),
                      ('connect2', 'conv_marker'),
                      ('connect2', 'convected'),
                      ('connect2', 'diffused')]
        assert str(ds) == "out.e-s002"
        assert ds.dimensionality == 3
        assert ds.current_time == 2.0
        assert ds.field_list == field_list

    #-----
    # test_gold
    #-----
    def test_gold(self):
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
        ds = utils.data_dir_load(gold)
        field_list = [('all', 'forced'), ('connect1', 'forced')]
        assert str(ds) == "gold.e"
        assert ds.field_list == field_list

    #-----
    # test_displacement_fields
    #-----
    def test_displacement_fields(self):
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
        displacement_dicts =[{'connect2': (5.0, [0.0, 0.0, 0.0])},
                             {'connect1': (1.0, [1.0, 2.0, 3.0]),
                              'connect2': (0.0, [0.0, 0.0, 0.0])}]
        ga_hd = b''
        for disp in displacement_dicts:
            ds = utils.data_dir_load(big_data, displacements=disp)
            for mesh in ds.index.meshes:
                def array_func():
                    return mesh.connectivity_coords
                ga_hd += self.generic_array_test(ds, array_func, 12)
        hashes = {}
        hashes['generic_array'] = utils.generate_hash(ga_hd)
        utils.handle_hashes(self.save_dir, 'exodusii-disp-fields', hashes, self.answer_store)

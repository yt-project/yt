"""
Title: test_tipsy.py
Purpose: Tipsy tests using the AGORA dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

from yt.testing import \
    assert_equal, \
    requires_file
from yt.frontends.tipsy.api import TipsyDataset
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = 'TipsyGalaxy/galaxy.00300'


# Globals
_fields = (("deposit", "all_density"),
           ("deposit", "all_count"),
           ("deposit", "DarkMatter_density"),
)
tg_fields = OrderedDict(
    [
        (('gas', 'density'), None),
        (('gas', 'temperature'), None),
        (('gas', 'temperature'), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (('gas', 'Fe_fraction'), None),
        (('Stars', 'Metals'), None),
    ]
)


#============================================
#                 TestTipsy
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestTipsy(fw.AnswerTest):
    """
    Container for tipsy frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_pkdgrav
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(pkdgrav, file_check = True)
    def test_pkdgrav(self, ds_pkdgrav):
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
        ds = ds_pkdgrav
        ppv_hd = b''
        fv_hd = b''
        dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
        dd = ds.all_data()
        assert_equal(dd["Coordinates"].shape, (26847360, 3))
        tot = sum(dd[ptype,"Coordinates"].shape[0]
                  for ptype in ds.particle_types if ptype != "all")
        assert_equal(tot, 26847360)
        for dobj_name in dso:
            for field in _fields:
                for axis in [0, 1, 2]:
                    for weight_field in [None]:
                        ppv_hd += self.pixelized_projection_values_test(
                            ds, axis, field, weight_field,
                            dobj_name)
                fv_hd += self.field_values_test(ds, field, dobj_name)
            dobj = utils.create_obj(ds, dobj_name)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert_equal(s1, s2)
        hashes = {'pixelized_projection_values' : utils.generate_hash(ppv_hd),
            'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'tipsy-test-pkdgrav', hashes, self.answer_store) 

    #-----
    # test_gasoline_dmonly
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(gasoline_dmonly, file_check = True)
    def test_gasoline_dmonly(self, ds_gasoline_dmonly):
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
        ds = ds_gasoline_dmonly
        ppv_hd = b''
        fv_hd = b''
        dso = [ None, ("sphere", ("c", (0.3, 'unitary')))]
        dd = ds.all_data()
        assert_equal(dd["Coordinates"].shape, (10550576, 3))
        tot = sum(dd[ptype,"Coordinates"].shape[0]
                  for ptype in ds.particle_types if ptype != "all")
        assert_equal(tot, 10550576)
        for dobj_name in dso:
            for field in _fields:
                for axis in [0, 1, 2]:
                    for weight_field in [None]:
                        ppv_hd += self.pixelized_projection_values_test(
                            ds, axis, field, weight_field,
                            dobj_name)
                fv_hd += self.field_values_test(ds, field, dobj_name)
            dobj = utils.create_obj(ds, dobj_name)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert_equal(s1, s2)
        hashes = {'pixelized_projection_values' : utils.generate_hash(ppv_hd),
            'field_values' : utils.generate_hash(fv_hd)}
        utils.handle_hashes(self.save_dir, 'tipsy-test-gasoline-dmonly', hashes, self.answer_store) 

    #-----
    # test_tipsy_galaxy
    #-----
    @utils.requires_ds(tipsy_gal)
    def test_tipsy_galaxy(self):
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
        ds = utils.data_dir_load(tipsy_gal)
        hashes = self.sph_answer(ds, 'galaxy.00300', 315372, tg_fields)
        utils.handle_hashes(self.save_dir, 'tipsy-test-galaxy', hashes, self.answer_store)

    #-----
    # test_TipsyDataset
    #-----
    @requires_file(gasoline_dmonly)
    @requires_file(pkdgrav)
    def test_TipsyDataset(self, ds_gasoline_dmonly, ds_pkdgrav):
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
        assert isinstance(ds_pkdgrav, TipsyDataset)
        assert isinstance(ds_gasoline_dmonly, TipsyDataset)

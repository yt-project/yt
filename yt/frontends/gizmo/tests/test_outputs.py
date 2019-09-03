"""
Title: test_gizmo.py
Purpose: Gizmo frontend tests
Notes:
    Copyright (c) 2015, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

import yt
from yt.testing import \
    requires_file
from yt.frontends.gizmo.api import GizmoDataset
from yt.frontends.gizmo.fields import metal_elements
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
g64 = "gizmo_64/output/snap_N64L16_135.hdf5"
gmhd = "gizmo_mhd_mwdisk/gizmo_mhd_mwdisk.hdf5"


# Globals
# This maps from field names to weight field names to use for projections
fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), ('gas', 'density')),
        (("gas", "metallicity"), ('gas', 'density')),
        (("gas", "O_metallicity"), ('gas', 'density')),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
    ]
)


#============================================
#                TestGizmo
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestGizmo(fw.AnswerTest):
    """
    Container for gizmo frontend answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """
    #-----
    # test_gizmo_64
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(g64)
    def test_gizmo_64(self):
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
        ds = yt.load(g64)
        assert isinstance(ds, GizmoDataset)
        hashes = self.sph_answer(ds, 'snap_N64L16_135', 524288, fields)
        utils.handle_hashes(self.save_dir, 'test-gizmo-64', hashes, self.answer_store)

    #-----
    # test_gizmo_mhd
    #-----
    @requires_file(gmhd)
    def test_gizmo_mhd(self, ds_gmhd):
        """
        Magnetic fields should be loaded correctly when they are present.

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
        ds = ds_gmhd
        ad = ds.all_data()
        ptype = 'PartType0'
        # Test vector magnetic field
        fmag = 'particle_magnetic_field'
        f = ad[ptype, fmag]
        assert str(f.units) == 'code_magnetic'
        assert f.shape == (409013, 3)
        # Test component magnetic fields
        for axis in 'xyz':
            f = ad[ptype, fmag + '_' + axis]
            assert str(f.units) == 'code_magnetic'
            assert f.shape == (409013,)

    #-----
    # test_gas_particle_fields
    #-----
    @requires_file(gmhd)
    def test_gas_particle_fields(self, ds_gmhd):
        """
        Test fields set up in GizmoFieldInfo.setup_gas_particle_fields.

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
        ds = ds_gmhd
        ptype = "PartType0"
        derived_fields = []
        # Add species fields
        for species in ["H", "H_p0", "H_p1"]:
            for suffix in ["density", "fraction", "mass", "number_density"]:
                derived_fields += ["%s_%s" % (species, suffix)]
        for species in metal_elements:
            derived_fields += ["%s_nuclei_mass_density" % species]
        # Add magnetic fields
        derived_fields += ["particle_magnetic_field_%s" % axis for axis in "xyz"]
        # Check
        for field in derived_fields:
            assert (ptype, field) in ds.derived_field_list
        ptype = "gas"
        derived_fields = []
        for species in ["H", "H_p0", "H_p1"]:
            for suffix in ["density", "number_density"]:
                derived_fields += ["%s_%s" % (species, suffix)]
        for species in metal_elements:
            for suffix in ["nuclei_mass_density", "metallicity"]:
                derived_fields += ["%s_%s" % (species, suffix)]
        for field in derived_fields:
            assert (ptype, field) in ds.derived_field_list

    #-----
    # test_star_particle_fields
    #-----
    @requires_file(gmhd)
    def test_star_particle_fields(self, ds_gmhd):
        """
        Test fields set up in GizmoFieldInfo.setup_star_particle_fields.

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
        ds = ds_gmhd
        ptype = "PartType4"
        derived_fields =[
            "creation_time",
            "age"
        ]
        for field in derived_fields:
            assert (ptype, field) in ds.derived_field_list

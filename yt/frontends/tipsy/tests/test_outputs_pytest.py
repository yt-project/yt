"""
Title: test_tipsy.py
Purpose: Tipsy tests using the AGORA dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.frontends.tipsy.api import TipsyDataset
from yt.testing import ParticleSelectionComparison, requires_file
from yt.utilities.answer_testing.answer_tests import nbody_answer, sph_answer
from yt.utilities.answer_testing.utils import data_dir_load, requires_ds

# Test data
pkdgrav = "halo1e11_run1.00400/halo1e11_run1.00400"
gasoline_dmonly = "agora_1e11.00400/agora_1e11.00400"
tipsy_gal = "TipsyGalaxy/galaxy.00300"


pkdgrav_cosmology_parameters = dict(
    current_redshift=0.0, omega_lambda=0.728, omega_matter=0.272, hubble_constant=0.702
)

pkdgrav_kwargs = dict(
    field_dtypes={"Coordinates": "d"},
    cosmology_parameters=pkdgrav_cosmology_parameters,
    unit_base={"length": (60.0, "Mpccm/h")},
)


@pytest.mark.answer_test
class TestTipsy:
    answer_file = None

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @requires_ds(pkdgrav, file_check=True)
    def test_pkdgrav(self, f, a, d, w):
        ds = data_dir_load(pkdgrav, TipsyDataset, (), kwargs=pkdgrav_kwargs)
        nb = nbody_answer(ds, "halo1e11_run1.00400", 26847360, f, w, d, a)
        self.hashes.update({"nbody_answer": nb})

    @pytest.mark.big_data
    @requires_ds(pkdgrav, file_check=True)
    def test_pkdgrav_psc(self):
        ds = data_dir_load(pkdgrav, TipsyDataset, (), kwargs=pkdgrav_kwargs)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.big_data
    @pytest.mark.usefixtures("hashing")
    @requires_ds(gasoline_dmonly, file_check=True)
    def test_gasoline_dmonly(self, f, a, d, w):
        cosmology_parameters = dict(
            current_redshift=0.0,
            omega_lambda=0.728,
            omega_matter=0.272,
            hubble_constant=0.702,
        )
        kwargs = dict(
            cosmology_parameters=cosmology_parameters,
            unit_base={"length": (60.0, "Mpccm/h")},
        )
        ds = data_dir_load(gasoline_dmonly, TipsyDataset, (), kwargs)
        nb = nbody_answer(ds, "agora_1e11.00400", 10550576, f, w, d, a)
        self.hashes.update({"nbody_answer": nb})

    @pytest.mark.big_data
    @requires_ds(gasoline_dmonly, file_check=True)
    def test_gasoline_dmonly_psc(self):
        cosmology_parameters = dict(
            current_redshift=0.0,
            omega_lambda=0.728,
            omega_matter=0.272,
            hubble_constant=0.702,
        )
        kwargs = dict(
            cosmology_parameters=cosmology_parameters,
            unit_base={"length": (60.0, "Mpccm/h")},
        )
        ds = data_dir_load(gasoline_dmonly, TipsyDataset, (), kwargs)
        psc = ParticleSelectionComparison(ds)
        psc.run_defaults()

    @pytest.mark.usefixtures("hashing")
    @requires_ds(tipsy_gal)
    def test_tipsy_galaxy_sph(self, f, w, d, a):
        ds = data_dir_load(
            tipsy_gal,
            kwargs={"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
        )
        # These tests should be re-enabled.  But the holdup is that the region
        # selector does not offset by domain_left_edge, and we have inelegant
        # selection using bboxes.
        # psc = ParticleSelectionComparison(ds)
        # psc.run_defaults()
        sph = sph_answer(ds, "galaxy.00300", 315372, f, w, d, a)
        self.hashes.update({"sph_answer": sph})

    @pytest.mark.usefixtures("hashing")
    @requires_ds(tipsy_gal)
    def test_tipsy_galaxy_nbody(self, f, w, d, a):
        ds = data_dir_load(
            tipsy_gal,
            kwargs={"bounding_box": [[-2000, 2000], [-2000, 2000], [-2000, 2000]]},
        )
        nb = nbody_answer(ds, "galaxy.00300", 315372, f, w, d, a)
        self.hashes.update({"nbody_answer": nb})

    @requires_file(gasoline_dmonly)
    @requires_file(pkdgrav)
    def test_TipsyDataset(self):
        assert isinstance(data_dir_load(pkdgrav, kwargs=pkdgrav_kwargs), TipsyDataset)
        assert isinstance(data_dir_load(gasoline_dmonly), TipsyDataset)

    @requires_file(tipsy_gal)
    def test_tipsy_index(self):
        ds = data_dir_load(tipsy_gal)
        sl = ds.slice("z", 0.0)
        assert sl["gas", "density"].shape[0] != 0

"""
Title: test_open_pmd.py
Purpose: openPMD frontend tests
Notes:
    Copyright (c) 2016, Fabian Koller (HZDR).
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from itertools import product

import numpy as np
import pytest

from yt.convenience import load
from yt.frontends.open_pmd.data_structures import OpenPMDDataset
from yt.testing import (
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
    requires_file,
)

# Test data
twoD = "example-2d/hdf5/data00000100.h5"
threeD = "example-3d/hdf5/data00000100.h5"
noFields = "no_fields/data00000400.h5"
noParticles = "no_particles/data00000400.h5"
groupBased = "singleParticle/simData.h5"


particle_fields = [
    "particle_charge",
    "particle_mass",
    "particle_momentum_x",
    "particle_momentum_y",
    "particle_momentum_z",
    "particle_positionCoarse_x",
    "particle_positionCoarse_y",
    "particle_positionCoarse_z",
    "particle_positionOffset_x",
    "particle_positionOffset_y",
    "particle_positionOffset_z",
    "particle_weighting",
]


@pytest.mark.answer_test
class TestOpenPMD:
    answer_file = None

    @pytest.mark.parametrize("ds", [threeD], indirect=True)
    def test_3d_out(self, ds):
        particle_types = ["all", "io", "nbody"]
        field_list = list(product(particle_types, particle_fields))
        field_list += list(product(("openPMD",), ("E_x", "E_y", "E_z", "rho")))
        domain_dimensions = [26, 26, 201] * np.ones_like(ds.domain_dimensions)
        domain_width = [2.08e-05, 2.08e-05, 2.01e-05] * np.ones_like(
            ds.domain_left_edge
        )
        assert isinstance(ds, OpenPMDDataset)
        assert_equal(ds.dimensionality, 3)
        assert_equal(ds.particle_types_raw, ("io",))
        assert "all" in ds.particle_unions
        assert_array_equal(ds.field_list, field_list)
        assert_array_equal(ds.domain_dimensions, domain_dimensions)
        assert_almost_equal(
            ds.current_time, 3.28471214521e-14 * np.ones_like(ds.current_time)
        )
        assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)

    @pytest.mark.parametrize("ds", [twoD], indirect=True)
    def test_2d_out(self, ds):
        particle_types = ("Hydrogen1+", "all", "electrons", "nbody")
        field_list = list(product(particle_types, particle_fields))
        field_list += [
            ("openPMD", "B_x"),
            ("openPMD", "B_y"),
            ("openPMD", "B_z"),
            ("openPMD", "E_x"),
            ("openPMD", "E_y"),
            ("openPMD", "E_z"),
            ("openPMD", "J_x"),
            ("openPMD", "J_y"),
            ("openPMD", "J_z"),
            ("openPMD", "rho"),
        ]
        domain_dimensions = [51, 201, 1] * np.ones_like(ds.domain_dimensions)
        domain_width = [3.06e-05, 2.01e-05, 1e0] * np.ones_like(ds.domain_left_edge)
        assert isinstance(ds, OpenPMDDataset)
        assert_equal(ds.dimensionality, 2)
        assert_equal(ds.particle_types_raw, ("Hydrogen1+", "electrons"))
        assert "all" in ds.particle_unions
        assert_array_equal(ds.field_list, field_list)
        assert_array_equal(ds.domain_dimensions, domain_dimensions)
        assert_almost_equal(
            ds.current_time, 3.29025596712e-14 * np.ones_like(ds.current_time)
        )
        assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)

    @pytest.mark.parametrize("ds", [noFields], indirect=True)
    def test_no_fields_out(self, ds):
        particle_types = ("all", "io", "nbody")
        no_fields_pfields = sorted(particle_fields + ["particle_id"])
        field_list = list(product(particle_types, no_fields_pfields))
        domain_dimensions = [1, 1, 1] * np.ones_like(ds.domain_dimensions)
        domain_width = [1, 1, 1] * np.ones_like(ds.domain_left_edge)
        assert isinstance(ds, OpenPMDDataset)
        assert_equal(ds.dimensionality, 3)
        assert_equal(ds.particle_types_raw, ("io",))
        assert "all" in ds.particle_unions
        assert_array_equal(ds.field_list, field_list)
        assert_array_equal(ds.domain_dimensions, domain_dimensions)
        assert_almost_equal(
            ds.current_time, 1.3161023868481013e-13 * np.ones_like(ds.current_time)
        )
        assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)

    @pytest.mark.parametrize("ds", [noParticles], indirect=True)
    def test_no_particles_out(self, ds):
        field_list = [
            ("openPMD", "E_x"),
            ("openPMD", "E_y"),
            ("openPMD", "E_z"),
            ("openPMD", "rho"),
        ]
        domain_dimensions = [51, 201, 1] * np.ones_like(ds.domain_dimensions)
        domain_width = [3.06e-05, 2.01e-05, 1e0] * np.ones_like(ds.domain_left_edge)
        assert isinstance(ds, OpenPMDDataset)
        assert_equal(ds.dimensionality, 2)
        assert_equal(ds.particle_types_raw, ("io",))
        assert "all" not in ds.particle_unions
        assert_array_equal(ds.field_list, field_list)
        assert_array_equal(ds.domain_dimensions, domain_dimensions)
        assert_almost_equal(
            ds.current_time, 1.3161023868481013e-13 * np.ones_like(ds.current_time)
        )
        assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)

    @requires_file(groupBased)
    def test_groupBased_out(self):
        dss = load(groupBased)
        particle_types = ("all", "io", "nbody")
        field_list = list(product(particle_types, particle_fields))
        field_list += [
            ("openPMD", "J_x"),
            ("openPMD", "J_y"),
            ("openPMD", "J_z"),
            ("openPMD", "e-chargeDensity"),
        ]
        domain_dimensions = [32, 64, 64] * np.ones_like(dss[0].domain_dimensions)
        domain_width = [0.0002752, 0.0005504, 0.0005504] * np.ones_like(
            dss[0].domain_left_edge
        )
        assert_equal(len(dss), 101)
        # Test only every 20th ds in the series
        for i in range(0, len(dss), 20):
            ds = dss[i]
            assert_equal(str(ds), "simData.h5")
            assert_equal(ds.dimensionality, 3)
            assert_equal(ds.particle_types_raw, ("io",))
            assert_array_equal(ds.field_list, field_list)
            assert_array_equal(ds.domain_dimensions, domain_dimensions)
            assert ds.current_time >= np.zeros_like(ds.current_time)
            assert ds.current_time <= 1.6499999999999998e-12 * np.ones_like(
                ds.current_time
            )
            assert_almost_equal(
                ds.domain_right_edge - ds.domain_left_edge, domain_width
            )

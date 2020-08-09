from itertools import product

import numpy as np

from yt.frontends.open_pmd.data_structures import OpenPMDDataset
from yt.loaders import load
from yt.testing import (
    assert_almost_equal,
    assert_array_equal,
    assert_equal,
    requires_file,
)
from yt.utilities.answer_testing.framework import data_dir_load

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


@requires_file(threeD)
def test_3d_out():
    ds = data_dir_load(threeD)
    particle_types = ["all", "io", "nbody"]
    field_list = list(product(particle_types, particle_fields))
    field_list += list(product(("openPMD",), ("E_x", "E_y", "E_z", "rho")))
    domain_dimensions = [26, 26, 201] * np.ones_like(ds.domain_dimensions)
    domain_width = [2.08e-05, 2.08e-05, 2.01e-05] * np.ones_like(ds.domain_left_edge)

    assert isinstance(ds, OpenPMDDataset)
    assert_equal(str(ds), "data00000100.h5")
    assert_equal(ds.dimensionality, 3)
    assert_equal(ds.particle_types_raw, ("io",))
    assert "all" in ds.particle_unions
    assert_array_equal(ds.field_list, field_list)
    assert_array_equal(ds.domain_dimensions, domain_dimensions)
    assert_almost_equal(
        ds.current_time, 3.28471214521e-14 * np.ones_like(ds.current_time)
    )
    assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)


@requires_file(twoD)
def test_2d_out():
    ds = data_dir_load(twoD)
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
    assert_equal(str(ds), "data00000100.h5")
    assert_equal(ds.dimensionality, 2)
    assert_equal(ds.particle_types_raw, ("Hydrogen1+", "electrons"))
    assert "all" in ds.particle_unions
    assert_array_equal(ds.field_list, field_list)
    assert_array_equal(ds.domain_dimensions, domain_dimensions)
    assert_almost_equal(
        ds.current_time, 3.29025596712e-14 * np.ones_like(ds.current_time)
    )
    assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)


@requires_file(noFields)
def test_no_fields_out():
    ds = data_dir_load(noFields)
    particle_types = ("all", "io", "nbody")
    no_fields_pfields = sorted(particle_fields + ["particle_id"])
    field_list = list(product(particle_types, no_fields_pfields))
    domain_dimensions = [1, 1, 1] * np.ones_like(ds.domain_dimensions)
    domain_width = [1, 1, 1] * np.ones_like(ds.domain_left_edge)

    assert isinstance(ds, OpenPMDDataset)
    assert_equal(str(ds), "data00000400.h5")
    assert_equal(ds.dimensionality, 3)
    assert_equal(ds.particle_types_raw, ("io",))
    assert "all" in ds.particle_unions
    assert_array_equal(ds.field_list, field_list)
    assert_array_equal(ds.domain_dimensions, domain_dimensions)
    assert_almost_equal(
        ds.current_time, 1.3161023868481013e-13 * np.ones_like(ds.current_time)
    )
    assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)


@requires_file(noParticles)
def test_no_particles_out():
    ds = data_dir_load(noParticles)
    field_list = [
        ("openPMD", "E_x"),
        ("openPMD", "E_y"),
        ("openPMD", "E_z"),
        ("openPMD", "rho"),
    ]
    domain_dimensions = [51, 201, 1] * np.ones_like(ds.domain_dimensions)
    domain_width = [3.06e-05, 2.01e-05, 1e0] * np.ones_like(ds.domain_left_edge)

    assert isinstance(ds, OpenPMDDataset)
    assert_equal(str(ds), "data00000400.h5")
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
def test_groupBased_out():
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
    for i in range(0, len(dss), 20):  # Test only every 20th ds out of the series
        ds = dss[i]
        assert_equal(str(ds), "simData.h5")
        assert_equal(ds.dimensionality, 3)
        assert_equal(ds.particle_types_raw, ("io",))
        assert_array_equal(ds.field_list, field_list)
        assert_array_equal(ds.domain_dimensions, domain_dimensions)
        assert ds.current_time >= np.zeros_like(ds.current_time)
        assert ds.current_time <= 1.6499999999999998e-12 * np.ones_like(ds.current_time)
        assert_almost_equal(ds.domain_right_edge - ds.domain_left_edge, domain_width)

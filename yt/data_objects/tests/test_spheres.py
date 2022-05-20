import numpy as np
from numpy.testing import assert_array_equal

from yt.data_objects.profiles import create_profile
from yt.testing import (
    assert_equal,
    assert_raises,
    fake_random_ds,
    periodicity_cases,
    requires_module,
)
from yt.utilities.exceptions import YTException, YTFieldNotFound


def setup():
    from yt.config import ytcfg

    ytcfg["yt", "internals", "within_testing"] = True


_fields_to_compare = (
    "spherical_radius",
    "cylindrical_radius",
    "spherical_theta",
    "cylindrical_theta",
    "spherical_phi",
    "cylindrical_z",
)


def test_domain_sphere():
    # Now we test that we can get different radial velocities based on field
    # parameters.

    # Get the first sphere
    ds = fake_random_ds(
        16,
        fields=("density", "velocity_x", "velocity_y", "velocity_z"),
        units=("g/cm**3", "cm/s", "cm/s", "cm/s"),
    )
    sp0 = ds.sphere(ds.domain_center, 0.25)

    # Compute the bulk velocity from the cells in this sphere
    bulk_vel = sp0.quantities.bulk_velocity()

    # Get the second sphere
    sp1 = ds.sphere(ds.domain_center, 0.25)

    # Set the bulk velocity field parameter
    sp1.set_field_parameter("bulk_velocity", bulk_vel)

    assert_equal(
        np.any(sp0[("gas", "radial_velocity")] == sp1[("gas", "radial_velocity")]),
        False,
    )

    # Radial profile without correction
    # Note we set n_bins = 8 here.

    rp0 = create_profile(
        sp0,
        "radius",
        "radial_velocity",
        units={"radius": "kpc"},
        logs={"radius": False},
        n_bins=8,
    )

    # Radial profile with correction for bulk velocity

    rp1 = create_profile(
        sp1,
        "radius",
        "radial_velocity",
        units={"radius": "kpc"},
        logs={"radius": False},
        n_bins=8,
    )

    assert_equal(rp0.x_bins, rp1.x_bins)
    assert_equal(rp0.used, rp1.used)
    assert_equal(rp0.used.sum() > rp0.used.size / 2.0, True)
    assert_equal(
        np.any(rp0["radial_velocity"][rp0.used] == rp1["radial_velocity"][rp1.used]),
        False,
    )

    ref_sp = ds.sphere("c", 0.25)
    for f in _fields_to_compare:
        ref_sp[f].sort()
    for center in periodicity_cases(ds):
        sp = ds.sphere(center, 0.25)
        for f in _fields_to_compare:
            sp[f].sort()
            assert_equal(sp[f], ref_sp[f])


def test_sphere_center():
    ds = fake_random_ds(
        16,
        nprocs=8,
        fields=("density", "temperature", "velocity_x"),
        units=("g/cm**3", "K", "cm/s"),
    )

    # Test if we obtain same center in different ways
    sp1 = ds.sphere("max", (0.25, "unitary"))
    sp2 = ds.sphere("max_density", (0.25, "unitary"))
    assert_array_equal(sp1.center, sp2.center)

    sp1 = ds.sphere("min", (0.25, "unitary"))
    sp2 = ds.sphere("min_density", (0.25, "unitary"))
    assert_array_equal(sp1.center, sp2.center)


def test_center_error():
    ds = fake_random_ds(16, nprocs=16)

    with assert_raises(YTFieldNotFound):
        ds.sphere("min_non_existing_field_name", (0.25, "unitary"))

    with assert_raises(YTFieldNotFound):
        ds.sphere("max_non_existing_field_name", (0.25, "unitary"))


@requires_module("miniball")
def test_minimal_sphere():
    ds = fake_random_ds(16, nprocs=8, particles=100)

    pos = ds.r[("all", "particle_position")]
    sp1 = ds.minimal_sphere(pos)

    N0 = len(pos)

    # Check all particles have been found
    N1 = len(sp1[("all", "particle_ones")])
    assert_equal(N0, N1)

    # Check that any smaller sphere is missing some particles
    sp2 = ds.sphere(sp1.center, sp1.radius * 0.9)
    N2 = len(sp2[("all", "particle_ones")])
    assert N2 < N0


@requires_module("miniball")
def test_minimal_sphere_bad_inputs():
    ds = fake_random_ds(16, nprocs=8, particles=100)
    pos = ds.r[("all", "particle_position")]

    ## Check number of points >= 2
    # -> should fail
    assert_raises(YTException, ds.minimal_sphere, pos[:1, :])

    # -> should not fail
    ds.minimal_sphere(pos[:2, :])

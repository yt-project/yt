import numpy as np

from yt.data_objects.profiles import create_profile
from yt.testing import \
    fake_random_ds, \
    assert_equal, \
    periodicity_cases

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

_fields_to_compare = ("spherical_r", "cylindrical_r",
                      "spherical_theta", "cylindrical_theta",
                      "spherical_phi", "cylindrical_z")

def test_domain_sphere():
    # Now we test that we can get different radial velocities based on field
    # parameters.

    # Get the first sphere
    ds = fake_random_ds(16, fields = ("density",
      "velocity_x", "velocity_y", "velocity_z"))
    sp0 = ds.sphere(ds.domain_center, 0.25)

    # Compute the bulk velocity from the cells in this sphere
    bulk_vel = sp0.quantities.bulk_velocity()

    # Get the second sphere
    sp1 = ds.sphere(ds.domain_center, 0.25)

    # Set the bulk velocity field parameter
    sp1.set_field_parameter("bulk_velocity", bulk_vel)

    yield assert_equal, np.any(sp0["radial_velocity"] ==
                               sp1["radial_velocity"]), False

    # Radial profile without correction
    # Note we set n_bins = 8 here.

    rp0 = create_profile(sp0, 'radius', 'radial_velocity',
                         units = {'radius': 'kpc'},
                         logs = {'radius': False},
                         n_bins = 8)

    # Radial profile with correction for bulk velocity

    rp1 = create_profile(sp1, 'radius', 'radial_velocity',
                         units = {'radius': 'kpc'},
                         logs = {'radius': False},
                         n_bins = 8)

    yield assert_equal, rp0.x_bins, rp1.x_bins
    yield assert_equal, rp0.used, rp1.used
    yield assert_equal, rp0.used.sum() > rp0.used.size/2.0, True
    yield assert_equal, np.any(rp0["radial_velocity"][rp0.used] ==
                               rp1["radial_velocity"][rp1.used]), \
                               False

    ref_sp = ds.sphere("c", 0.25)
    for f in _fields_to_compare:
        ref_sp[f].sort()
    for center in periodicity_cases(ds):
        sp = ds.sphere(center, 0.25)
        for f in _fields_to_compare:
            sp[f].sort()
            yield assert_equal, sp[f], ref_sp[f]

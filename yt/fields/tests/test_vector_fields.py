import numpy as np

from yt.testing import assert_allclose_units, fake_random_ds, requires_file
from yt.units import cm, s
from yt.utilities.answer_testing.framework import data_dir_load


def random_unit_vector(prng):
    v = prng.random_sample(3)
    while (v == 0).all():
        v = prng.random_sample(3)
    return v / np.sqrt((v ** 2).sum())


def random_velocity_vector(prng):
    return 2e5 * prng.random_sample(3) - 1e5


def compare_vector_conversions(data_source):
    prng = np.random.RandomState(8675309)
    normals = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] + [
        random_unit_vector(prng) for i in range(2)
    ]
    bulk_velocities = [random_velocity_vector(prng) for i in range(2)]

    for bv in bulk_velocities:
        bulk_velocity = bv * cm / s
        data_source.set_field_parameter("bulk_velocity", bulk_velocity)
        data_source.clear_data()

        vmag = data_source["velocity_magnitude"]
        vrad = data_source["velocity_spherical_radius"]

        for normal in normals:
            data_source.set_field_parameter("normal", normal)
            data_source.clear_data()

            assert_allclose_units(vrad, data_source["velocity_spherical_radius"])

            vmag_new = data_source["velocity_magnitude"]
            assert_allclose_units(vmag, vmag_new)

            vmag_cart = np.sqrt(
                (data_source["velocity_x"] - bulk_velocity[0]) ** 2
                + (data_source["velocity_y"] - bulk_velocity[1]) ** 2
                + (data_source["velocity_z"] - bulk_velocity[2]) ** 2
            )
            assert_allclose_units(vmag, vmag_cart)

            vmag_cyl = np.sqrt(
                data_source["velocity_cylindrical_radius"] ** 2
                + data_source["velocity_cylindrical_theta"] ** 2
                + data_source["velocity_cylindrical_z"] ** 2
            )
            assert_allclose_units(vmag, vmag_cyl)

            vmag_sph = np.sqrt(
                data_source["velocity_spherical_radius"] ** 2
                + data_source["velocity_spherical_theta"] ** 2
                + data_source["velocity_spherical_phi"] ** 2
            )
            assert_allclose_units(vmag, vmag_sph)

            for i, d in enumerate("xyz"):
                assert_allclose_units(
                    data_source[f"velocity_{d}"] - bulk_velocity[i],
                    data_source[f"relative_velocity_{d}"],
                )

        for i, ax in enumerate("xyz"):
            data_source.set_field_parameter("axis", i)
            data_source.clear_data()
            assert_allclose_units(
                data_source["velocity_los"], data_source[f"relative_velocity_{ax}"]
            )

        for i, ax in enumerate("xyz"):
            prj = data_source.ds.proj("velocity_los", i, weight_field="density")
            assert_allclose_units(prj["velocity_los"], prj[f"velocity_{ax}"])

        data_source.clear_data()
        ax = [0.1, 0.2, -0.3]
        data_source.set_field_parameter("axis", ax)
        ax /= np.sqrt(np.dot(ax, ax))
        vlos = data_source["relative_velocity_x"] * ax[0]
        vlos += data_source["relative_velocity_y"] * ax[1]
        vlos += data_source["relative_velocity_z"] * ax[2]
        assert_allclose_units(data_source["velocity_los"], vlos)


def test_vector_component_conversions_fake():
    ds = fake_random_ds(16)
    ad = ds.all_data()
    compare_vector_conversions(ad)


g30 = "IsolatedGalaxy/galaxy0030/galaxy0030"


@requires_file(g30)
def test_vector_component_conversions_real():
    ds = data_dir_load(g30)
    sp = ds.sphere(ds.domain_center, (10, "kpc"))
    compare_vector_conversions(sp)

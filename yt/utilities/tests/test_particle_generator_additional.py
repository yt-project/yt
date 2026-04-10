import numpy as np
from numpy.testing import assert_allclose, assert_array_equal, assert_raises

from yt.loaders import load_uniform_grid
from yt.utilities.particle_generator import (
    FromListParticleGenerator,
    LatticeParticleGenerator,
    ParticleGenerator,
)


def _make_uniform_ds(nprocs=1):
    density = np.ones((4, 4, 4))
    temperature = 4.0 * np.ones((4, 4, 4))
    fields = {
        "density": (density, "code_mass/code_length**3"),
        "temperature": (temperature, "K"),
    }
    return load_uniform_grid(fields, density.shape, 1.0, nprocs=nprocs)


def _field_list(*extra_fields):
    return [
        ("io", "particle_position_x"),
        ("io", "particle_position_y"),
        ("io", "particle_position_z"),
        *extra_fields,
    ]


def test_particle_generator_base_accessors_and_setup_fields():
    ds = _make_uniform_ds()
    particles = ParticleGenerator(ds, 2, _field_list(("io", "particle_mass")))

    assert len(particles) == 2
    assert particles.has_key(("io", "particle_mass"))
    assert particles.keys() == _field_list(("io", "particle_mass")) + [
        ("io", "particle_index")
    ]

    x = np.array([0.2, 0.4])
    y = np.array([0.3, 0.5])
    z = np.array([0.1, 0.2])
    mass = np.array([1.5, 2.5])
    setup_fields = {
        ("io", "particle_position_x"): np.array([9.0, 8.0]),
        ("io", "particle_mass"): mass,
    }

    particles._setup_particles(x, y, z, setup_fields=setup_fields)

    assert_allclose(particles["io", "particle_position_x"], x)
    assert_allclose(particles["io", "particle_position_y"], y)
    assert_allclose(particles["io", "particle_position_z"], z)
    assert_allclose(particles["io", "particle_mass"], mass)

    updated_mass = np.array([3.5, 4.5])
    particles["io", "particle_mass"] = updated_mass
    assert_allclose(particles["io", "particle_mass"], updated_mass)


def test_particle_generator_requires_all_position_fields():
    ds = _make_uniform_ds()

    assert_raises(
        KeyError,
        ParticleGenerator,
        ds,
        1,
        [("io", "particle_position_x"), ("io", "particle_position_y")],
    )


def test_from_list_particle_generator_string_keys_and_domain_errors():
    ds = _make_uniform_ds()

    valid_data = {
        "particle_position_x": np.array([0.2, 0.4]),
        "particle_position_y": np.array([0.3, 0.5]),
        "particle_position_z": np.array([0.1, 0.2]),
        "particle_mass": np.array([5.0, 6.0]),
    }
    particles = FromListParticleGenerator(ds, 2, valid_data.copy())
    assert_allclose(particles["io", "particle_mass"], [5.0, 6.0])
    assert_allclose(particles["io", "particle_position_x"], [0.2, 0.4])

    outside_domain_data = {
        "particle_position_x": np.array([1.0]),
        "particle_position_y": np.array([0.25]),
        "particle_position_z": np.array([0.25]),
    }
    assert_raises(
        ValueError,
        FromListParticleGenerator,
        ds,
        1,
        outside_domain_data.copy(),
    )

    missing_positions_data = {("io", "particle_mass"): np.array([7.0])}
    assert_raises(
        UnboundLocalError,
        FromListParticleGenerator,
        ds,
        1,
        missing_positions_data.copy(),
    )


def test_lattice_particle_generator_rejects_domain_edge_bounds():
    ds = _make_uniform_ds()

    assert_raises(
        ValueError,
        LatticeParticleGenerator,
        ds,
        np.array([2, 2, 2]),
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.5, 0.5]),
        _field_list(),
    )


def test_particle_generator_multigrid_empty_grid_paths_and_density_mapping():
    ds = _make_uniform_ds(nprocs=2)
    particles = ParticleGenerator(
        ds,
        2,
        _field_list(("io", "particle_gas_density")),
    )

    particles._setup_particles(
        np.array([0.25, 0.75]),
        np.array([0.25, 0.75]),
        np.array([0.10, 0.40]),
    )

    assert_array_equal(particles.NumberOfParticles, [2, 0])
    assert_array_equal(particles.ParticleGridIndices, [0, 2, 2])

    particles.map_grid_fields_to_particles(
        {("gas", "density"): ("io", "particle_gas_density")}
    )
    assert_allclose(particles["io", "particle_gas_density"], [1.0, 1.0])

    particles.apply_to_stream(overwrite=True)

    grid0, grid1 = ds.index.grids
    assert grid0["all", "particle_position_x"].size == 2
    assert grid1["all", "particle_position_x"].size == 0
    assert_allclose(grid0["all", "particle_gas_density"], [1.0, 1.0])

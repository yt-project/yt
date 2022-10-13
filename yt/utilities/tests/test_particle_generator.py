import numpy as np

from yt.loaders import load_uniform_grid
from yt.testing import assert_almost_equal, assert_equal
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.utilities.particle_generator import (
    FromListParticleGenerator,
    LatticeParticleGenerator,
    WithDensityParticleGenerator,
)


def test_particle_generator():
    # First generate our dataset
    domain_dims = (32, 32, 32)
    dens = np.zeros(domain_dims) + 0.1
    temp = 4.0 * np.ones(domain_dims)
    fields = {"density": (dens, "code_mass/code_length**3"), "temperature": (temp, "K")}
    ds = load_uniform_grid(fields, domain_dims, 1.0)

    # Now generate particles from density

    field_list = [
        ("io", "particle_position_x"),
        ("io", "particle_position_y"),
        ("io", "particle_position_z"),
        ("io", "particle_index"),
        ("io", "particle_gas_density"),
    ]
    num_particles = 10000
    field_dict = {("gas", "density"): ("io", "particle_gas_density")}
    sphere = ds.sphere(ds.domain_center, 0.45)

    particles1 = WithDensityParticleGenerator(ds, sphere, num_particles, field_list)
    particles1.assign_indices()
    particles1.map_grid_fields_to_particles(field_dict)

    # Test to make sure we ended up with the right number of particles per grid
    particles1.apply_to_stream()
    particles_per_grid1 = [grid.NumberOfParticles for grid in ds.index.grids]
    assert_equal(particles_per_grid1, particles1.NumberOfParticles)
    particles_per_grid1 = [
        len(grid[("all", "particle_position_x")]) for grid in ds.index.grids
    ]
    assert_equal(particles_per_grid1, particles1.NumberOfParticles)

    tags = uconcatenate([grid[("all", "particle_index")] for grid in ds.index.grids])
    assert np.unique(tags).size == num_particles

    del tags

    # Set up a lattice of particles
    pdims = np.array([32, 32, 32])

    def new_indices():
        # We just add new indices onto the existing ones
        return np.arange(np.product(pdims)) + num_particles

    le = np.array([0.25, 0.25, 0.25])
    re = np.array([0.75, 0.75, 0.75])

    particles2 = LatticeParticleGenerator(ds, pdims, le, re, field_list)
    particles2.assign_indices(function=new_indices)
    particles2.map_grid_fields_to_particles(field_dict)

    # Test lattice positions
    xpos = np.unique(particles2["io", "particle_position_x"])
    ypos = np.unique(particles2["io", "particle_position_y"])
    zpos = np.unique(particles2["io", "particle_position_z"])

    xpred = np.linspace(le[0], re[0], num=pdims[0], endpoint=True)
    ypred = np.linspace(le[1], re[1], num=pdims[1], endpoint=True)
    zpred = np.linspace(le[2], re[2], num=pdims[2], endpoint=True)

    assert_almost_equal(xpos, xpred)
    assert_almost_equal(ypos, ypred)
    assert_almost_equal(zpos, zpred)

    del xpos, ypos, zpos
    del xpred, ypred, zpred

    # Test the number of particles again
    particles2.apply_to_stream()
    particles_per_grid2 = [grid.NumberOfParticles for grid in ds.index.grids]
    assert_equal(
        particles_per_grid2, particles1.NumberOfParticles + particles2.NumberOfParticles
    )

    [grid.field_data.clear() for grid in ds.index.grids]
    particles_per_grid2 = [
        len(grid[("all", "particle_position_x")]) for grid in ds.index.grids
    ]
    assert_equal(
        particles_per_grid2, particles1.NumberOfParticles + particles2.NumberOfParticles
    )

    # Test the uniqueness of tags
    tags = np.concatenate([grid[("all", "particle_index")] for grid in ds.index.grids])
    tags.sort()
    assert_equal(tags, np.arange(np.product(pdims) + num_particles))

    del tags

    # Now dump all of these particle fields out into a dict
    pdata = {}
    dd = ds.all_data()
    for field in field_list:
        pdata[field] = dd[field]

    # Test the "from-list" generator and particle field overwrite
    num_particles3 = num_particles + np.product(pdims)
    particles3 = FromListParticleGenerator(ds, num_particles3, pdata)
    particles3.apply_to_stream(overwrite=True)

    # Test the number of particles again
    particles_per_grid3 = [grid.NumberOfParticles for grid in ds.index.grids]
    assert_equal(
        particles_per_grid3, particles1.NumberOfParticles + particles2.NumberOfParticles
    )
    particles_per_grid2 = [
        len(grid[("all", "particle_position_z")]) for grid in ds.index.grids
    ]
    assert_equal(
        particles_per_grid3, particles1.NumberOfParticles + particles2.NumberOfParticles
    )
    assert_equal(particles_per_grid2, particles_per_grid3)

    # Test adding in particles with a different particle type

    num_star_particles = 20000
    pdata2 = {
        ("star", "particle_position_x"): np.random.uniform(size=num_star_particles),
        ("star", "particle_position_y"): np.random.uniform(size=num_star_particles),
        ("star", "particle_position_z"): np.random.uniform(size=num_star_particles),
    }

    particles4 = FromListParticleGenerator(ds, num_star_particles, pdata2, ptype="star")
    particles4.apply_to_stream()

    dd = ds.all_data()
    assert dd["star", "particle_position_x"].size == num_star_particles
    assert dd["io", "particle_position_x"].size == num_particles3
    assert dd["all", "particle_position_x"].size == num_star_particles + num_particles3

    del pdata
    del pdata2
    del ds
    del particles1
    del particles2
    del particles4
    del fields
    del dens
    del temp

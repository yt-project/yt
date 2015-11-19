import numpy as np
from yt.testing import \
    assert_almost_equal, \
    assert_equal
from yt.utilities.particle_generator import \
    WithDensityParticleGenerator, \
    LatticeParticleGenerator, \
    FromListParticleGenerator
from yt.frontends.stream.api import load_uniform_grid, refine_amr
import yt.utilities.initial_conditions as ic
import yt.utilities.flagging_methods as fm
from yt.units.yt_array import uconcatenate

def setup() :
    pass

def test_particle_generator():
    # First generate our dataset
    domain_dims = (128, 128, 128)
    dens = np.zeros(domain_dims) + 0.1
    temp = 4.*np.ones(domain_dims)
    fields = {"density": (dens, 'code_mass/code_length**3'),
              "temperature": (temp, 'K')}
    ug = load_uniform_grid(fields, domain_dims, 1.0)
    fo = [ic.BetaModelSphere(1.0,0.1,0.5,[0.5,0.5,0.5],{"density":(10.0)})]
    rc = [fm.flagging_method_registry["overdensity"](4.0)]
    ds = refine_amr(ug, rc, fo, 3)

    # Now generate particles from density

    field_list = [("io", "particle_position_x"),
                  ("io", "particle_position_y"),
                  ("io", "particle_position_z"),
                  ("io", "particle_index"),
                  ("io", "particle_gas_density")]
    num_particles = 1000000
    field_dict = {("gas", "density"): ("io", "particle_gas_density")}
    sphere = ds.sphere(ds.domain_center, 0.45)

    particles1 = WithDensityParticleGenerator(ds, sphere, num_particles, field_list)
    particles1.assign_indices()
    particles1.map_grid_fields_to_particles(field_dict)
    
    # Test to make sure we ended up with the right number of particles per grid
    particles1.apply_to_stream()
    particles_per_grid1 = [grid.NumberOfParticles for grid in ds.index.grids]
    yield assert_equal, particles_per_grid1, particles1.NumberOfParticles
    particles_per_grid1 = [len(grid["particle_position_x"]) for grid in ds.index.grids]
    yield assert_equal, particles_per_grid1, particles1.NumberOfParticles

    tags = uconcatenate([grid["particle_index"] for grid in ds.index.grids])
    assert(np.unique(tags).size == num_particles)
    # Set up a lattice of particles
    pdims = np.array([64,64,64])
    def new_indices() :
        # We just add new indices onto the existing ones
        return np.arange((np.product(pdims)))+num_particles
    le = np.array([0.25,0.25,0.25])
    re = np.array([0.75,0.75,0.75])
    new_field_list = field_list + [("io", "particle_gas_temperature")]
    new_field_dict = {("gas", "density"): ("io", "particle_gas_density"),
                      ("gas", "temperature"): ("io", "particle_gas_temperature")}

    particles2 = LatticeParticleGenerator(ds, pdims, le, re, new_field_list)
    particles2.assign_indices(function=new_indices)
    particles2.map_grid_fields_to_particles(new_field_dict)

    #Test lattice positions
    xpos = np.unique(particles2["io", "particle_position_x"])
    ypos = np.unique(particles2["io", "particle_position_y"])
    zpos = np.unique(particles2["io", "particle_position_z"])

    xpred = np.linspace(le[0],re[0],num=pdims[0],endpoint=True)
    ypred = np.linspace(le[1],re[1],num=pdims[1],endpoint=True)
    zpred = np.linspace(le[2],re[2],num=pdims[2],endpoint=True)

    assert_almost_equal( xpos, xpred)
    assert_almost_equal( ypos, ypred)
    assert_almost_equal( zpos, zpred)

    #Test the number of particles again
    particles2.apply_to_stream()
    particles_per_grid2 = [grid.NumberOfParticles for grid in ds.index.grids]
    yield assert_equal, particles_per_grid2, particles1.NumberOfParticles+particles2.NumberOfParticles

    [grid.field_data.clear() for grid in ds.index.grids]
    particles_per_grid2 = [len(grid["particle_position_x"]) for grid in ds.index.grids]
    yield assert_equal, particles_per_grid2, particles1.NumberOfParticles+particles2.NumberOfParticles

    #Test the uniqueness of tags
    tags = np.concatenate([grid["particle_index"] for grid in ds.index.grids])
    tags.sort()
    yield assert_equal, tags, np.arange((np.product(pdims)+num_particles))

    # Test that the old particles have zero for the new field
    old_particle_temps = [grid["particle_gas_temperature"][:particles_per_grid1[i]]
                          for i, grid in enumerate(ds.index.grids)]
    test_zeros = [np.zeros((particles_per_grid1[i])) 
                  for i, grid in enumerate(ds.index.grids)]
    yield assert_equal, old_particle_temps, test_zeros

    #Now dump all of these particle fields out into a dict
    pdata = {}
    dd = ds.all_data()
    for field in new_field_list :
        pdata[field] = dd[field]

    #Test the "from-list" generator and particle field clobber
    particles3 = FromListParticleGenerator(ds, num_particles+np.product(pdims), pdata)
    particles3.apply_to_stream(clobber=True)
    
    #Test the number of particles again
    particles_per_grid3 = [grid.NumberOfParticles for grid in ds.index.grids]
    yield assert_equal, particles_per_grid3, particles1.NumberOfParticles+particles2.NumberOfParticles
    particles_per_grid2 = [len(grid["particle_position_z"]) for grid in ds.index.grids]
    yield assert_equal, particles_per_grid3, particles1.NumberOfParticles+particles2.NumberOfParticles

if __name__=="__main__":
    for n, i in enumerate(test_particle_generator()):
        i[0](*i[1:])

import numpy as np
from yt.mods import *
from yt.testing import *
from yt.frontends.stream.api import load_uniform_grid, refine_amr, load_amr_grids
import yt.utilities.initial_conditions as ic
import yt.utilities.flagging_methods as fm

def setup() :
    pass

# Field information

def test_stream_particles() :
    num_particles = 100000
    domain_dims = (64, 64, 64)
    dens = np.random.random(domain_dims) 
    x = np.random.uniform(size=num_particles)
    y = np.random.uniform(size=num_particles)
    z = np.random.uniform(size=num_particles)
    m = np.ones((num_particles))

    # Field operators and cell flagging methods

    fo = []
    fo.append(ic.TopHatSphere(0.1, [0.2,0.3,0.4],{"Density": 2.0}))
    fo.append(ic.TopHatSphere(0.05, [0.7,0.4,0.75],{"Density": 20.0}))
    rc = [fm.flagging_method_registry["overdensity"](1.0)]
    
    # Check that all of this runs ok without particles
    
    ug0 = load_uniform_grid({"Density": dens}, domain_dims, 1.0)
    ug0 = load_uniform_grid({"Density": dens}, domain_dims, 1.0, nprocs=8)
    amr0 = refine_amr(ug0, rc, fo, 3)

    grid_data = []
    
    for grid in amr0.h.grids :
        
        data = dict(left_edge = grid.LeftEdge,
                    right_edge = grid.RightEdge,
                    level = grid.Level,
                    dimensions = grid.ActiveDimensions,
                    number_of_particles = grid.NumberOfParticles)
    
        for field in amr0.h.field_list :
            
            data[field] = grid[field]
            
        grid_data.append(data)

    amr0 = load_amr_grids(grid_data, domain_dims, 1.0)
                        
    # Now add particles

    fields1 = {"Density": dens,
               "particle_position_x": x,
               "particle_position_y": y,
               "particle_position_z": z,
               "particle_mass": m,
               "number_of_particles": num_particles}

    fields2 = fields1.copy()

    ug1 = load_uniform_grid(fields1, domain_dims, 1.0)
    ug2 = load_uniform_grid(fields2, domain_dims, 1.0, nprocs=8)

    # Check to make sure the number of particles is the same

    number_of_particles1 = np.sum([grid.NumberOfParticles for grid in ug1.h.grids])
    number_of_particles2 = np.sum([grid.NumberOfParticles for grid in ug2.h.grids])
    
    yield assert_equal, number_of_particles1, num_particles
    yield assert_equal, number_of_particles1, number_of_particles2

    # Check to make sure the fields have been defined correctly
    
    assert ug1._get_field_info("all", "particle_position_x").particle_type
    assert ug1._get_field_info("all", "particle_position_y").particle_type
    assert ug1._get_field_info("all", "particle_position_z").particle_type
    assert ug1._get_field_info("all", "particle_mass").particle_type
    assert not ug1._get_field_info("gas", "Density").particle_type

    assert ug2._get_field_info("all", "particle_position_x").particle_type
    assert ug2._get_field_info("all", "particle_position_y").particle_type
    assert ug2._get_field_info("all", "particle_position_z").particle_type
    assert ug2._get_field_info("all", "particle_mass").particle_type
    assert not ug2._get_field_info("gas", "Density").particle_type
    
    # Now refine this

    amr1 = refine_amr(ug1, rc, fo, 3)
    
    grid_data = []
    
    for grid in amr1.h.grids :
        
        data = dict(left_edge = grid.LeftEdge,
                    right_edge = grid.RightEdge,
                    level = grid.Level,
                    dimensions = grid.ActiveDimensions,
                    number_of_particles = grid.NumberOfParticles)

        for field in amr1.h.field_list :

            data[field] = grid[field]
            
        grid_data.append(data)
    
    amr2 = load_amr_grids(grid_data, domain_dims, 1.0)

    # Check everything again

    number_of_particles1 = [grid.NumberOfParticles for grid in amr1.h.grids]
    number_of_particles2 = [grid.NumberOfParticles for grid in amr2.h.grids]
    
    yield assert_equal, np.sum(number_of_particles1), num_particles
    yield assert_equal, number_of_particles1, number_of_particles2
    
    assert amr1._get_field_info("all", "particle_position_x").particle_type
    assert amr1._get_field_info("all", "particle_position_y").particle_type
    assert amr1._get_field_info("all", "particle_position_z").particle_type
    assert amr1._get_field_info("all", "particle_mass").particle_type
    assert not amr1._get_field_info("gas", "Density").particle_type
    
    assert amr2._get_field_info("all", "particle_position_x").particle_type
    assert amr2._get_field_info("all", "particle_position_y").particle_type
    assert amr2._get_field_info("all", "particle_position_z").particle_type
    assert amr2._get_field_info("all", "particle_mass").particle_type
    assert not amr2._get_field_info("gas", "Density").particle_type

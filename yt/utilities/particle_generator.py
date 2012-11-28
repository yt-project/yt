import numpy as np
import h5py
from yt.utilities.lib import CICSample_3
from yt.funcs import *

class ParticleGenerator(object) :

    default_fields = ["particle_position_x",
                      "particle_position_y",
                      "particle_position_z",
                      "particle_index"]

    def __init__(self, pf, num_particles, field_list) :
        
        self.pf = pf
        self.num_particles = num_particles
        self.field_list = field_list
            
        try :
            self.posx_index = self.field_list.index(self.default_fields[0])
            self.posy_index = self.field_list.index(self.default_fields[1])
            self.posz_index = self.field_list.index(self.default_fields[2])
            self.index_index = self.field_list.index(self.default_fields[3])
        except :
            raise KeyError("Field list must contain the following fields: " +
                           "\'particle_position_x\', \'particle_position_y\'" +
                           ", \'particle_position_z\', \'particle_index\' ")

        self.num_grids = self.pf.h.num_grids
        self.NumberOfParticles = np.zeros((self.num_grids), dtype='int64')
        self.ParticleIndices = np.zeros(self.num_grids + 1, dtype='int64')
        
        self.num_fields = len(self.field_list)
        
        self.particles = np.zeros((self.num_particles, self.num_fields),
                                  dtype='float64')

    def has_key(self, key) :

        return (key in self.field_list)
            
    def keys(self) :

        return self.field_list.keys()
    
    def __getitem__(self, key) :

        """
        Get the field associated with key.
        """

        return self.particles[:,self.field_list.index(key)]
    
    def __setitem__(self, key, val) :
        
        """
        Sets a field to be some other value.
        """

        self.particles[:,self.field_list.index(key)] = val[:]
                
    def __len__(self) :

        """
        The number of particles
        """
        
        return self.num_particles

    def get_from_grid(self, grid) :

        ind = grid.id-grid._id_offset
        
        start = self.ParticleIndices[ind]
        end = self.ParticleIndices[ind+1]

        return dict([(field, self.particles[start:end,self.field_list.index(field)])
                     for field in self.field_list])
    
    def setup_particles(self,x,y,z) :

        particle_grids, particle_grid_inds = self.pf.h.find_points(x,y,z)

        idxs = np.argsort(particle_grid_inds)

        self.particles[:,self.posx_index] = x[idxs]
        self.particles[:,self.posy_index] = y[idxs]
        self.particles[:,self.posz_index] = z[idxs]

        self.NumberOfParticles = np.bincount(particle_grid_inds,
                                             minlength=self.num_grids)

        if self.num_grids > 1 :
            np.add.accumulate(self.NumberOfParticles.squeeze(),
                              out=self.ParticleIndices[1:])
        else :
            self.ParticleIndices[1] = self.NumberOfParticles.squeeze()

        return idxs
    
    def assign_indices(self, function=None, **kwargs) :

        if function is None :
            self.particles[:,self.index_index] = np.arange((self.num_particles))
        else :
            self.particles[:,self.index_index] = function(**kwargs)
            
    def map_grid_fields_to_particles(self, mapping_dict) :
        
        pbar = get_pbar("Mapping fields to particles", self.num_grids)
        
        for i, grid in enumerate(self.pf.h.grids) :

            pbar.update(i)
            
            if self.NumberOfParticles[i] >  0:

                start = self.ParticleIndices[i]
                end = self.ParticleIndices[i+1]

                # Note we add one ghost zone!
                
                cube = grid.retrieve_ghost_zones(1, mapping_dict.keys())

                le = np.array(grid.LeftEdge).astype(np.float64)
                dims = np.array(grid.ActiveDimensions).astype(np.int32)
                
                for gfield, pfield in mapping_dict.items() :

                    field_index = self.field_list.index(pfield)
                    
                    CICSample_3(self.particles[start:end,self.posx_index],
                                self.particles[start:end,self.posy_index],
                                self.particles[start:end,self.posz_index],
                                self.particles[start:end,field_index],
                                np.int64(self.NumberOfParticles[i]),
                                cube[gfield], le, dims,
                                np.float64(grid['dx']))

        pbar.finish()

    def apply_to_stream(self) :

        grid_data = []
        
        for i,g in enumerate(self.pf.h.grids) :

            data = {}
                        
            data["number_of_particles"] = self.NumberOfParticles[i] + \
                                          g.NumberOfParticles

            grid_particles = self.get_from_grid(g)
            
            for field in self.field_list :

                if data["number_of_particles"] > 0 :
                
                    if g.NumberOfParticles > 0 :
                        data[field] = np.concatenate(g[field],
                                                     grid_particles[field])
                    else :
                        data[field] = grid_particles[field]
                        
                else :

                    data[field] = np.array([], dtype='float64')

            grid_data.append(data)
            
        self.pf.h.update_data(grid_data)

class FromListParticleGenerator(ParticleGenerator) :

    def __init__(self, pf, num_particles, data) :

        field_list = data.keys()

        x = data["particle_position_x"]
        y = data["particle_position_y"]
        z = data["particle_position_z"]

        xcond = np.logical_or(x < pf.domain_left_edge[0],
                              x >= pf.domain_right_edge[0])
        ycond = np.logical_or(y < pf.domain_left_edge[1],
                              y >= pf.domain_right_edge[1])
        zcond = np.logical_or(z < pf.domain_left_edge[2],
                              z >= pf.domain_right_edge[2])

        cond = np.logical_or(xcond, ycond)
        cond = np.logical_or(zcond, cond)

        if np.any(cond) :
            raise ValueError("Some particles are outside of the domain!!!")
            
        ParticleGenerator.__init__(self, pf, num_particles, field_list)

        idxs = self.setup_particles(x,y,z)
        
        for field in field_list :
            
            self.particles[:,self.field_list.index(field)] = data[field][idxs]

class LatticeParticleGenerator(ParticleGenerator) :

    def __init__(self, pf, particles_dims, particles_left_edge,
                 particles_right_edge, field_list) :

        num_x = particles_dims[0]
        num_y = particles_dims[1]
        num_z = particles_dims[2]
        
        xmin = particles_left_edge[0]
        ymin = particles_left_edge[1]
        zmin = particles_left_edge[2]
        
        xmax = particles_right_edge[0]
        ymax = particles_right_edge[1]
        zmax = particles_right_edge[2]
                                
        xcond = (xmin < pf.domain_left_edge[0]) or \
                (xmax >= pf.domain_right_edge[0])
        ycond = (ymin < pf.domain_left_edge[1]) or \
                (ymax >= pf.domain_right_edge[1])
        zcond = (zmin < pf.domain_left_edge[2]) or \
                (zmax >= pf.domain_right_edge[2])
                
        cond = xcond or ycond or zcond

        if cond :
            raise ValueError("Proposed bounds for particles are outside domain!!!")
        
        ParticleGenerator.__init__(self, pf, num_x*num_y*num_z, field_list)

        dx = (xmax-xmin)/(num_x-1)
        dy = (ymax-ymin)/(num_y-1)
        dz = (zmax-zmin)/(num_z-1)
        
        inds = np.indices((num_x,num_y,num_z))
        xpos = inds[0]*dx + xmin
        ypos = inds[1]*dy + ymin
        zpos = inds[2]*dz + zmin
        
        self.setup_particles(xpos.flat, ypos.flat, zpos.flat)
        
class WithDensityParticleGenerator(ParticleGenerator) :

    def __init__(self, pf, data_source, num_particles, field_list,
                 density_field="Density") :

        ParticleGenerator.__init__(self, pf, num_particles, field_list)

        num_cells = len(data_source["x"].flat)
        
        max_density = data_source[density_field].max()

        num_particles_left = num_particles

        all_x = []
        all_y = []
        all_z = []
        
        pbar = get_pbar("Generating Particles", num_particles)

        tot_num_accepted = int(0)
        
        while num_particles_left > 0:

            rho = np.random.uniform(high=1.01*max_density,
                                    size=num_particles_left)

            idxs = np.random.random_integers(low=0, high=num_cells-1,
                                             size=num_particles_left)

            rho_true = data_source[density_field].flat[idxs]
            accept = rho <= rho_true
            num_accepted = accept.sum()

            accepted_idxs = idxs[accept]
            
            xpos = data_source["x"].flat[accepted_idxs] + \
                   np.random.uniform(low=-0.5, high=0.5, size=num_accepted) * \
                   data_source["dx"].flat[accepted_idxs]
            ypos = data_source["y"].flat[accepted_idxs] + \
                   np.random.uniform(low=-0.5, high=0.5, size=num_accepted) * \
                   data_source["dy"].flat[accepted_idxs]
            zpos = data_source["z"].flat[accepted_idxs] + \
                   np.random.uniform(low=-0.5, high=0.5, size=num_accepted) * \
                   data_source["dz"].flat[accepted_idxs]

            all_x.append(xpos)
            all_y.append(ypos)
            all_z.append(zpos)

            num_particles_left -= num_accepted
            tot_num_accepted += num_accepted
            
            pbar.update(tot_num_accepted)

        pbar.finish()

        x = np.concatenate(all_x)
        y = np.concatenate(all_y)
        z = np.concatenate(all_z)

        self.setup_particles(x,y,z)
        

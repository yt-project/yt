"""
Find a progenitor line by reverse traversing a timeseries
and finding the max density around the previous timestep

Author: Christopher Moody <chrisemoody@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

#Take the first snapshot, find max
#Define new position for 2nd snap from old position + bulk velocity
#Draw sphere around new position
#Find new max density
#Or particle approach:
#Take last snapshot
#Find highest particle density
#Keep indices of n=1000 most bound particles
#Then go to N-1 snapshot, find new positions of particles
#Find new max particle density 
#Find new set of n=1000 particles around new center

class FindMaxProgenitor:
    def __init__(self,ts):
        self.ts = ts
    
    def find_max_field(self,field,radius=0.01)):
        """
        Find the maximum of the given field, and iterate backwards through 
        snapshots finding the maxima within a small radius of the previous 
        snapshot's center.
        """
        centers = []
        v,c = ts[-1].h.find_max(field)
        for pf in ts[::-1]:
            sph = pf.h.sphere(c,radius)
            vn,c = sph.find_max(field)
            centers.append(c)
        return centers

    def find_max_particle(self,radius=0.01,nparticles=1000):
        """
        Find the particle at the maximum density and iterate backwards through
        snapshots, finding the localtion of the maximum density of the 
        previously closest nparticles.
        """
        indices = None
        dd = ts[-1].h.all_data()
        domain_dimensions = ts[-1].domain_dimensions
        sim_unit_to_cm = ts[-1]['cm']
        v,c = dd.quantities["ParticleDensityCenter"]()
        centers = [c,]
        for pf in ts[::-1]:
            if indices:
                dd = pf.h.all_data()
                data = dict(number_of_particles=indices.shape[0])
                for ax in 'xyz':
                    data['particle_position_%s'%ax]=\
                            dd["particle_position_%s"%ax][indices]
                subselection = load_uniform_grid(data,domain_dimensions,
                                                 sim_unit_to_cm)
                ss = subselection.h.all_data()
                v,c = ss.quantities["ParticleDensityCenter"]()
            sph = pf.h.sphere(c,radius)
            rad = sph["particle_radius"]
            idx = sph["particle_index"]
            indices = idx[np.argsort(rad)[:nparticles]]
            centers.append(c)
        return centers
            

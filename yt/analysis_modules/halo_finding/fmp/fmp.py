"""
Find a progenitor line by reverse traversing a timeseries
and finding the max density around the previous timestep

Author: Christopher Moody <chrisemoody@gmail.com>
Affiliation: UC Santa Cruz
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
import numpy as np
from yt.frontends.stream.api import load_uniform_grid

class FindMaxProgenitor:
    def __init__(self,ts):
        self.ts = ts
    
    def find_max_field(self,field="Density",radius=0.01,use_bulk=True):
        """
        Find the maximum of the given field, and iterate backwards through 
        snapshots finding the maxima within a small radius of the previous 
        snapshot's center.
        """
        centers = []
        v,c = self.ts[-1].h.find_max(field)
        t_old = None
        for pf in self.ts[::-1]:
            t = pf.current_time
            if t_old:
                dt = t_old - t
                c += dt*bv
            sph = pf.h.sphere(c,radius)
            v,i,x,y,z=sph.quantities["MaxLocation"]("Density")
            c = np.array([x,y,z])
            centers.append(c)
            bv = sph.quantities["BulkVelocity"]()
            #bv is in cgs but center is in untary
            bv /= pf['cm']
            t_old = pf.current_time
        return centers

    def find_max_particle(self,initial_center=None,radius=0.01,nparticles=1000,
                          particle_type="all"):
        """
        Find the particle at the maximum density and iterate backwards through
        snapshots, finding the localtion of the maximum density of the 
        previously closest nparticles.
        """
        indices = None
        dd = self.ts[-1].h.all_data()
        domain_dimensions = self.ts[-1].domain_dimensions
        sim_unit_to_cm = self.ts[-1]['cm']
        c = initial_center
        if c is None:
            v,c = dd.quantities["ParticleDensityCenter"](particle_type=\
                                                         particle_type)
        centers = []
        for pf in self.ts[::-1]:
            if indices is not None:
                dd = pf.h.sphere(c,radius)
                data = dict(number_of_particles=indices.shape[0])
                index = dd[(particle_type,'particle_index')]
                inside = find_index_in_array(index,indices)
                for ax in 'xyz':
                    pos = dd[(particle_type,"particle_position_%s"%ax)][inside]
                    data[('all','particle_position_%s'%ax)]= pos
                mas = dd[(particle_type,"particle_mass")][inside]
                data[('all','particle_mass')]= mas
                subselection = load_uniform_grid(data,domain_dimensions,
                                                 sim_unit_to_cm)
                ss = subselection.h.all_data()
                v,c = ss.quantities["ParticleDensityCenter"]()
            sph = pf.h.sphere(c,radius)
            rad = sph["ParticleRadius"]
            idx = sph["particle_index"]
            indices = idx[np.argsort(rad)[:nparticles]]
            centers.append(c)
        return centers
            
def chunks(l, n):
    #http://stackoverflow.com/questions/312443/how-do-you-split-
    #a-list-into-evenly-sized-chunks-in-python
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def find_index_in_array(arr1,arr2,size=long(1e6)):
    #for element in arr2 find corresponding index in arr1
    #temporary size is arr1.shape x arr2.shape so chunk this out
    indices = np.array((),'i8')
    for chunk in chunks(arr1,size):
        idx = np.where(np.reshape(chunk,(chunk.shape[0],1))==
                       np.reshape(arr2,(1,arr2.shape[0])))[0]
        indices = np.concatenate((indices,idx)).astype('i8')
    return indices


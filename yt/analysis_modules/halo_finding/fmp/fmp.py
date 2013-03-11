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
from yt.funcs import *
from yt.frontends.stream.api import load_uniform_grid
from yt.utilities.parallel_tools import parallel_analysis_interface
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     parallel_objects
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, ProcessorPool, Communicator

class FindMaxProgenitor(ParallelAnalysisInterface):
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
        centers = {}
        earliest_c = c
        for pfs_chunk in chunks(self.ts[::-1]):
            pf = pfs_chunk[0] #first is the latest in time
            mylog.info("Finding central indices")
            sph = pf.h.sphere(earliest_c,radius)
            rad = sph["ParticleRadius"]
            idx = sph["particle_index"]
            indices = idx[np.argsort(rad)[:nparticles]]
            for sto, pf in parallel_objects(pfs_chunk,storage=centers):
                dd = pf.h.sphere(c,radius)
                data = dict(number_of_particles=indices.shape[0])
                index = dd[(particle_type,'particle_index')]
                inside = np.in1d(indices,index,assume_unique=True)
                mylog.info("Collecting particles %1.1e of %1.1e",inside.sum(),
                           nparticles)
                if inside.sum()==0:
                    mylog.warning("Found no matching indices in %s",str(pf)) 
                for ax in 'xyz':
                    pos = dd[(particle_type,"particle_position_%s"%ax)][inside]
                    data[('all','particle_position_%s'%ax)]= pos
                mas = dd[(particle_type,"particle_mass")][inside]
                data[('all','particle_mass')]= mas
                mylog.info("Finding center")
                subselection = load_uniform_grid(data,domain_dimensions,
                                                 sim_unit_to_cm)
                ss = subselection.h.all_data()
                v,c = ss.quantities["ParticleDensityCenter"]()
                sto.result_id = pf.parameters['aexpn']
                sto.result = c
            #last in the chunk is the earliest in time
            earliest_c = centers[pf_chunk[-1].parameters['aexpn']]
        return centers

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    #http://stackoverflow.com/questions/312443/
    #how-do-you-split-a-list-into-evenly-sized-chunks-in-python
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


from yt.data_objects.data_containers import YTFieldData
from yt.data_objects.time_series import TimeSeriesData
from yt.funcs import *

import numpy as na
import h5py

class ParticleTrajectoryCollection(object) :

    def __init__(self, filenames, indices, fields = None) :

        indices.sort() # Just in case the caller wasn't careful
        
        self.field_data = YTFieldData()
        self.pfs = TimeSeriesData.from_filenames(filenames)
        self.masks = []
        self.sorts = []
        self.indices = indices
        self.num_indices = len(indices)
        self.num_steps = len(filenames)
        self.times = []
        
        if fields is None : fields = ["particle_position_x",
                                      "particle_position_y",
                                      "particle_position_z"]

        for pf in self.pfs :
            dd = pf.h.all_data()
            newtags = dd["particle_index"].astype("int")
            if not na.all(na.in1d(indices, newtags)) :
                print "Not all requested particle ids contained in this file!"
                raise IndexError
            mask = na.in1d(newtags, indices, assume_unique=True)
            sorts = na.argsort(newtags[mask])
            self.masks.append(mask)            
            self.sorts.append(sorts)
            self.times.append(pf.current_time)

        self.times = na.array(self.times)
        
        self.particle_fields = [field for field in self.pfs[0].h.derived_field_list
                                if (field.startswith("particle") or
                                    field.startswith("Particle"))]

        for field in fields :

            self.get_data(field)
            
    def has_key(self, key) :

        return (key in self.field_data)
    
    def keys(self) :

        return self.field_data.keys()

    def __getitem__(self, key) :
            
        if key not in self.particle_fields :
            print "Not a valid particle field!"
            raise KeyError
        
        if not self.field_data.has_key(key) :

            self.get_data(key)

        return self.field_data[key]
    
    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.field_data[key] = val
                        
    def __delitem__(self, key) :

        del self.field_data[key]

    def __iter__(self) :
        
        for idx in xrange(self.num_indices) :
            traj = {}
            traj["particle_index"] = self.indices[idx]
            traj["particle_time"] = self.times
            for field in self.field_data.keys() :
                traj[field] = self[field][idx,:]
            yield traj
            
    def __len__(self) :

        return self.num_indices

    def add_fields(self, fields) :

        for field in fields :

            if not self.field_data.has_key(field):

                self.get_data(field)
                
    def get_data(self, field) :
            
        if not self.field_data.has_key(field):
            
            particles = na.empty((0))
                            
            for pf, mask, sort in zip(self.pfs, self.masks, self.sorts) :
            
                dd = pf.h.all_data()
                pfield = dd[field][mask]
                particles = na.append(particles, pfield[sort])
            
            self[field] = particles.reshape(self.num_steps,
                                            self.num_indices).transpose()

        return self.field_data[field]

    def trajectory_from_index(self, index) :
            
        mask = na.in1d(self.indices, (index,), assume_unique=True)

        if not na.any(mask) :
            print "The particle index %d is not in the list!" % (index)
            raise IndexError

        fields = [field for field in sorted(self.field_data.keys())]
                                
        traj = {}

        traj["particle_time"] = self.times
        traj["particle_index"] = index
        
        for field in fields :

            traj[field] = self[field][mask,:][0]

        return traj

    def write_out(self, filename_base) :

        fields = [field for field in sorted(self.field_data.keys())]

        num_fields = len(fields)

        first_str = "# particle_time\t" + "\t".join(fields)+"\n"
        
        template_str = "%g\t"*num_fields+"%g\n"
        
        for ix in xrange(self.num_indices) :

            outlines = [first_str]

            for it in xrange(self.num_steps) :
                outlines.append(template_str %
                                tuple([self.times[it]]+[self[field][ix,it] for field in fields]))
            
            fid = open(filename_base + "_%d.dat" % self.indices[ix], "w")
            fid.writelines(outlines)
            fid.close()
            del fid
            
    def write_out_h5(self, filename) :

        fid = h5py.File(filename, "w")

        fields = [field for field in sorted(self.field_data.keys())]
        
        fid.create_dataset("particle_indices", dtype=na.int32,
                           data=self.indices)
        fid.create_dataset("particle_time", data=self.times)
        
        for field in fields :

            fid.create_dataset("%s" % field, data=self[field])
                        
        fid.close()

"""
Author: John ZuHone <jzuhone@gmail.com>
Affiliation: NASA/GSFC
License:
  Copyright (C) 2012 John ZuHone All Rights Reserved.

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

from yt.data_objects.data_containers import YTFieldData
from yt.data_objects.time_series import TimeSeriesData
from yt.utilities.lib import sample_field_at_positions
from yt.convenience import load
from yt.funcs import *

import numpy as na
import h5py

class ParticleTrajectoryCollection(object) :

    r"""A collection of particle trajectories in time over a series of
    parameter files. 

    The ParticleTrajectoryCollection object contains a collection of
    particle trajectories for a specified set of particle indices. 
    
    Parameters
    ----------
    filenames : list of strings
        A time-sorted list of filenames to construct the TimeSeriesData
        object.
    indices : array_like
        An integer array of particle indices whose trajectories we
        want to track. If they are not sorted they will be sorted.
    fields : list of strings, optional
        A set of fields that is retrieved when the trajectory
        collection is instantiated.
        Default : None (will default to the fields 'particle_position_x',
        'particle_position_y', 'particle_position_z')

    Examples
    ________
    >>> from yt.mods import *
    >>> my_fns = glob.glob("orbit_hdf5_chk_00[0-9][0-9]")
    >>> my_fns.sort()
    >>> fields = ["particle_position_x", "particle_position_y",
    >>>           "particle_position_z", "particle_velocity_x",
    >>>           "particle_velocity_y", "particle_velocity_z"]
    >>> pf = load(my_fns[0])
    >>> init_sphere = pf.h.sphere(pf.domain_center, (.5, "unitary"))
    >>> indices = init_sphere["particle_index"].astype("int")
    >>> trajs = ParticleTrajectoryCollection(my_fns, indices, fields=fields)
    >>> for t in trajs :
    >>>     print t["particle_velocity_x"].max(), t["particle_velocity_x"].min()

    Notes
    -----
    As of this time only particle trajectories that are complete over the
    set of specified parameter files are supported. If any particle's history
    ends for some reason (e.g. leaving the simulation domain or being actively
    destroyed), the whole trajectory collection of which it is a set must end
    at or before the particle's last timestep. This is a limitation we hope to
    lift at some point in the future.     
    """
    def __init__(self, filenames, indices, fields = None) :

        indices.sort() # Just in case the caller wasn't careful
        
        self.field_data = YTFieldData()
        #self.pfs = TimeSeriesData.from_filenames(filenames)
        self.pfs = [load(fn) for fn in filenames]
        self.masks = []
        self.sorts = []
        self.indices = indices
        self.num_indices = len(indices)
        self.num_steps = len(filenames)
        self.times = []

        # Default fields 
        
        if fields is None : fields = []

        # Must ALWAYS have these fields
        
        fields = fields + ["particle_position_x",
                           "particle_position_y",
                           "particle_position_z"]

        """
        The following loops through the parameter files
        and performs two tasks. The first is to isolate
        the particles with the correct indices, and the
        second is to create a sorted list of these particles.
        We also make a list of the current time from each file. 
        Right now, the code assumes (and checks for) the
        particle indices existing in each file, a limitation I
        would like to lift at some point since some codes
        (e.g., FLASH) destroy particles leaving the domain.
        """
        
        for pf in self.pfs :
            dd = pf.h.all_data()
            newtags = dd["particle_index"].astype("int")
            #if not na.all(na.in1d(indices, newtags, assume_unique=True)) :
            #    print "Not all requested particle ids contained in this file!"
            #    raise IndexError
            #mask = na.in1d(newtags, indices, assume_unique=True)
            #sorts = na.argsort(newtags[mask])
            #self.masks.append(mask)            
            #self.sorts.append(sorts)
            #self.times.append(pf.current_time)

        self.times = na.array(self.times)

        # Set up the derived field list and the particle field list
        # so that if the requested field is a particle field, we'll
        # just copy the field over, but if the field is a grid field,
        # we will first copy the field over to the particle positions
        # and then return the field. 

        self.derived_field_list = self.pfs[0].h.derived_field_list
        self.particle_fields = [field for field in self.derived_field_list
                                if self.pfs[0].field_info[field].particle_type]

        # Now instantiate the requested fields 
        for field in fields :

            self._get_data(field)
            
    def has_key(self, key) :

        return (key in self.field_data)
    
    def keys(self) :

        return self.field_data.keys()

    def __getitem__(self, key) :
        """
        Get the field associated with key,
        checking to make sure it is a particle field.
        """

        if not self.field_data.has_key(key) :

            self._get_data(key)

        return self.field_data[key]
    
    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.field_data[key] = val
                        
    def __delitem__(self, key) :
        """
        Delete the field from the trajectory
        """
        del self.field_data[key]

    def __iter__(self) :

        """
        This iterates over the trajectories for
        the different particles, returning dicts
        of fields for each trajectory
        """
        for idx in xrange(self.num_indices) :
            traj = {}
            traj["particle_index"] = self.indices[idx]
            traj["particle_time"] = self.times
            for field in self.field_data.keys() :
                traj[field] = self[field][idx,:]
            yield traj
            
    def __len__(self) :

        """
        The number of individual trajectories
        """
        return self.num_indices

    def add_fields(self, fields) :

        """
        Add a list of fields to an existing trajectory

        Parameters
        ----------
        fields : list of strings
            A list of fields to be added to the current trajectory
            collection.

        Examples
        ________
        >>> from yt.mods import *
        >>> trajs = ParticleTrajectoryCollection(my_fns, indices)
        >>> trajs.add_fields(["particle_mass", "particle_gpot"])
        """
        
        for field in fields :

            if not self.field_data.has_key(field):

                self._get_data(field)
                
    def _get_data(self, field) :

        """
        Get a field to include in the trajectory collection.
        The trajectory collection itself is a dict of 2D numpy arrays,
        with shape (num_indices, num_steps)
        """
        
        if not self.field_data.has_key(field):
            
            particles = na.empty((0))

            step = int(0)
                
            for pf, mask, sort in zip(self.pfs, self.masks, self.sorts) :
                                    
                if field in self.particle_fields :

                    # This is easy... just get the particle fields

                    dd = pf.h.all_data()
                    pfield = dd[field][mask]
                    particles = na.append(particles, pfield[sort])

                else :

                    # This is hard... must loop over grids

                    pfield = na.zeros((self.num_indices))
                    x = self["particle_position_x"][:,step]
                    y = self["particle_position_y"][:,step]
                    z = self["particle_position_z"][:,step]

                    leaf_grids = [g for g in pf.h.grids if len(g.Children) == 0]
                        
                    for grid in leaf_grids :

                        pfield += sample_field_at_positions(grid[field],
                                                            grid.LeftEdge,
                                                            grid.RightEdge,
                                                            x, y, z)

                    particles = na.append(particles, pfield)

                step += 1
                
            self[field] = particles.reshape(self.num_steps,
                                            self.num_indices).transpose()

        return self.field_data[field]

    def trajectory_from_index(self, index) :

        """
        Retrieve a single trajectory corresponding to a specific particle
        index

        Parameters
        ----------
        index : int
            This defines which particle trajectory from the
            ParticleTrajectoryCollection object will be returned.

        Returns
        -------
        A dictionary corresponding to the particle's trajectory and the
        fields along that trajectory

        Examples
        --------
        >>> from yt.mods import *
        >>> import matplotlib.pylab as pl
        >>> trajs = ParticleTrajectoryCollection(my_fns, indices)
        >>> traj = trajs.trajectory_from_index(indices[0])
        >>> pl.plot(traj["particle_time"], traj["particle_position_x"], "-x")
        >>> pl.savefig("orbit")
        """
        
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

        """
        Write out particle trajectories to tab-separated ASCII files (one
        for each trajectory) with the field names in the file header. Each
        file is named with a basename and the index number.

        Parameters
        ----------
        filename_base : string
            The prefix for the outputted ASCII files.

        Examples
        --------
        >>> from yt.mods import *
        >>> trajs = ParticleTrajectoryCollection(my_fns, indices)
        >>> trajs.write_out("orbit_trajectory")       
        """
        
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

        """
        Write out all the particle trajectories to a single HDF5 file
        that contains the indices, the times, and the 2D array for each
        field individually

        Parameters
        ---------
        filename : string
            The output filename for the HDF5 file

        Examples
        --------
        >>> from yt.mods import *
        >>> trajs = ParticleTrajectoryCollection(my_fns, indices)
        >>> trajs.write_out_h5("orbit_trajectories")                
        """
        
        fid = h5py.File(filename, "w")

        fields = [field for field in sorted(self.field_data.keys())]
        
        fid.create_dataset("particle_indices", dtype=na.int32,
                           data=self.indices)
        fid.create_dataset("particle_time", data=self.times)
        
        for field in fields :

            fid.create_dataset("%s" % field, data=self[field])
                        
        fid.close()

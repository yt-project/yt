"""
FOF-output data handling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos.fof import *

class FOFList(object):
    def __init__(self, data_source, link=0.2,
                 dm_only = True):
        """
        Run fof on *data_source* with a given link length *link*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles. It's probably best to not run on star particles because
        friends-of-friends doesn't account for mass. 
        Returns an iterable collection of *FOFGroup* items.
        """
        self.data_source = data_source
        self.dm_only = dm_only
        self.link = link
        self._groups = []
        mylog.info("Initializing FOF")
        self.__obtain_particles()
        self.__run_fof()
        mylog.info("Parsing outputs")
        self.__parse_output()
        mylog.debug("Finished. (%s)", len(self))

    def __obtain_particles(self):
        if self.dm_only: ii = self.__get_dm_indices()
        else: ii = slice(None)
        self.particle_fields = {}
        for field in ["particle_position_%s" % ax for ax in 'xyz']:
            tot_part = self.data_source[field].size
            self.particle_fields[field] = self.data_source[field][ii]
        self._base_indices = na.arange(tot_part)[ii]

    def __run_fof(self):
        self.tags = \
            RunFOF(self.particle_fields["particle_position_x"],
                   self.particle_fields["particle_position_y"],
                   self.particle_fields["particle_position_z"],
                   self.link)
        self.particle_fields["tags"] = self.tags

    def __get_dm_indices(self):
        if 'creation_time' in self.data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on creation time")
            return (self.data_source["creation_time"] < 0)
        elif 'particle_type' in self.data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on particle type")
            return (self.data_source["particle_type"] == 1)
        else:
            mylog.warning("No particle_type, no creation_time, so not distinguishing.")
            return slice(None)

    def __parse_output(self):
        unique_ids = na.unique(self.tags)
        counts = na.bincount(self.tags+1)
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i+1]
            if i == -1:
                cp += counts[i+1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(FOFGroup(self, i, group_indices))
            px, py, pz = [self.particle_fields['particle_position_%s'%ax][group_indices]
                                            for ax in 'xyz']
            cp += counts[i+1]

    def __len__(self):
        return len(self._groups)
 
    def __iter__(self):
        return FOFIterator(self)

    def __getitem__(self, key):
        return self._groups[key]

    def write_out(self, filename="FOFAnalysis.out"):
        """
        Write out standard FOF information to *filename*.
        """
        if hasattr(filename, 'write'):
            f = filename
        else:
            f = open(filename,"w")
        f.write("\t".join(["   # Group","Mass","         # part",
                           "center-of-mass",
                           "x","y","z",
                           "vx","vy","vz","max_r","\n"]))
        for group in self:
            f.write("%10i\t" % group.id)
            f.write("%0.9e\t" % group.total_mass())
            f.write("%10i\t" % group.get_size())
            f.write("\t".join(["%0.9e" % v for v in group.center_of_mass()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.bulk_velocity()]))
            f.write("\t")
            f.write("%0.9e\t" % group.maximum_radius())
            f.write("\n")
            f.flush()
        f.close()

class FOFIterator(object):
    def __init__(self, fof):
        self.fof = fof
        self.index = -1

    def next(self):
        self.index += 1
        if self.index == len(self.fof): raise StopIteration
        return self.fof[self.index]

class FOFGroup(object):
    """
    A data source that returns particle information about the members of a
    FOF-identified halo.
    """
    __metaclass__ = ParallelDummy # This will proxy up our methods
    _distributed = False
    _processing = False
    _owner = 0
    indices = None
    dont_wrap = ["get_sphere", "write_particle_list"]
    extra_wrap = ["__getitem__"]

    def __init__(self, fof_output, id, indices = None):
        self.fof_output = fof_output
        self.id = id
        self.data = fof_output.data_source
        if indices is not None: self.indices = fof_output._base_indices[indices]
        # We assume that if indices = None, the instantiator has OTHER plans
        # for us -- i.e., setting it somehow else

    def center_of_mass(self):
        """
        Calculate and return the center of mass.
        """
        pm = self["ParticleMassMsun"]
        cx = self["particle_position_x"]
        cy = self["particle_position_y"]
        cz = self["particle_position_z"]
        c_vec = na.array([cx[0],cy[0],cz[0]]) - na.array([0.5,0.5,0.5])
        cx = cx - c_vec[0]
        cy = cy - c_vec[1]
        cz = cz - c_vec[2]
        com = na.array([v-na.floor(v) for v in [cx,cy,cz]])
        com = (pm * com).sum(axis=1)/pm.sum() + c_vec
        return com

    def total_mass(self):
        """
        Returns the total mass in solar masses of the halo.
        """
        return self["ParticleMassMsun"].sum()

    def bulk_velocity(self):
        """
        Returns the mass-weighted average velocity.
        """
        pm = self["ParticleMassMsun"]
        vx = (self["particle_velocity_x"] * pm).sum()
        vy = (self["particle_velocity_y"] * pm).sum()
        vz = (self["particle_velocity_z"] * pm).sum()
        return na.array([vx,vy,vz])/pm.sum()

    def maximum_radius(self):
        """
        Returns the maximum radius in the halo for all particles,
        either from the center_of_mass.
        """
        center = self.center_of_mass()
        rx = na.abs(self["particle_position_x"]-center[0])
        ry = na.abs(self["particle_position_y"]-center[1])
        rz = na.abs(self["particle_position_z"]-center[2])
        r = na.sqrt(na.minimum(rx, 1.0-rx)**2.0
                +   na.minimum(ry, 1.0-ry)**2.0
                +   na.minimum(rz, 1.0-rz)**2.0)
        return r.max()

    def __getitem__(self, key):
        return self.data[key][self.indices]

    def get_sphere(self):
        """
        Returns an EnzoSphere centered on 
        the *center_of_mass*, with the maximum radius of the halo.
        """
        center = self.center_of_mass()
        radius = self.maximum_radius()
        # A bit of a long-reach here...
        sphere = self.fof_output.data_source.hierarchy.sphere(
                        center, radius=radius)
        return sphere

    def get_size(self):
        return self.indices.size

    def write_particle_list(self, handle):
        self._processing = True
        gn = "Halo%08i" % (self.id)
        handle.createGroup("/", gn)
        for field in ["particle_position_%s" % ax for ax in 'xyz'] \
                   + ["particle_velocity_%s" % ax for ax in 'xyz'] \
                   + ["particle_index"]:
            handle.createArray("/%s" % gn, field, self[field])
        n = handle.getNode("/", gn)
        # set attributes on n
        self._processing = False

class FOFHaloFinder(FOFList, ParallelAnalysisInterface):
    def __init__(self, pf, link=0.2, dm_only=True, padding=0.2):
        self.pf = pf
        self.hierarchy = pf.h
        self.center = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0
        self.padding = 0.0 #* pf["unitary"] # This should be clevererer
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self.data_source = self._partition_hierarchy_3d(padding=self.padding)
        n_parts = self._mpi_allsum(self.data_source["particle_position_x"].size)
        print 'n_parts %d' % n_parts
        # get the average spacing between particles
        l = pf["DomainRightEdge"] - pf["DomainLeftEdge"]
        vol = l[0] * l[1] * l[2]
        avg_spacing = (float(vol) / n_parts)**(1./3.)
        print 'avg_spacing %f' % avg_spacing
        self.padding = padding
        padded, LE, RE, self.data_source = self._partition_hierarchy_3d(padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        self._reposition_particles((LE, RE))
        self.data_source.get_data(["particle_position_%s" % ax for ax in 'xyz'])
        # here is where the FOF halo finder is run
        super(FOFHaloFinder, self).__init__(self.data_source, link * avg_spacing, dm_only)
        self._parse_foflist()
        self._join_foflists()

    def _parse_foflist(self):
        groups, com, hi  = [], {}, 0
        LE, RE = self.bounds
        for halo in self._groups:
            this_com = halo.center_of_mass()
            # if the center of mass is in the box, keep it
            if na.all((this_com >= LE) & (this_com <= RE)):
                # Now we add the halo information to OURSELVES, taken from the
                # self.fof_list
                groups.append(FOFGroup(self, hi))
                groups[-1].indices = halo.indices
                self._claim_object(groups[-1])
                hi += 1
        del self._groups # explicit >> implicit
        self._groups = groups

    def _join_foflists(self):
        # First we get the total number of halos the entire collection
        # has identified
        # Note I have added a new method here to help us get information
        # about processors and ownership and so forth.
        # _mpi_info_dict returns a dict of {proc: whatever} where whatever is
        # what is fed in on each proc.
        mine, halo_info = self._mpi_info_dict(len(self))
        nhalos = sum(halo_info.values())
        # Figure out our offset
        my_first_id = sum([v for k,v in halo_info.items() if k < mine])
        # sort the list by the size of the groups
        # Now we add ghost halos and reassign all the IDs
        # Note: we already know which halos we own!
        after = my_first_id + len(self._groups)
        # One single fake halo, not owned, does the trick
        self._groups = [FOFGroup(self, i) for i in range(my_first_id)] + \
                       self._groups + \
                       [FOFGroup(self, i) for i in range(after, nhalos)]
        id = 0
        for proc in sorted(halo_info.keys()):
            for halo in self._groups[id:id+halo_info[proc]]:
                halo.id = id
                halo._distributed = self._distributed
                halo._owner = proc
                id += 1
        self._groups.sort(key = lambda h: -1 * h.get_size())
        for j,halo in enumerate(self._groups):
            halo.id = j
        
    def _reposition_particles(self, bounds):
        # This only does periodicity.  We do NOT want to deal with anything
        # else.  The only reason we even do periodicity is the 
        LE, RE = bounds
        dw = self.pf["DomainRightEdge"] - self.pf["DomainLeftEdge"]
        for i, ax in enumerate('xyz'):
            arr = self.data_source["particle_position_%s" % ax]
            arr[arr < LE[i]-self.padding] += dw[i]
            arr[arr > RE[i]+self.padding] -= dw[i]

    def write_out(self, filename):
        f = self._write_on_root(filename)
        FOFList.write_out(self, f)

    @parallel_blocking_call
    def write_particle_lists(self, prefix):
        fn = "%s.h5" % self._get_filename(prefix)
        f = tables.openFile(fn, "w")
        for halo in self._groups:
            if not self._is_mine(halo): continue
            halo.write_particle_list(f)

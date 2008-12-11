"""
HOP-output data handling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

from yt.lagos.hop import *

class HopList(object):
    def __init__(self, data_source, threshold=160.0,
                 dm_only = True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self.data_source = data_source
        self.dm_only = dm_only
        self.threshold = threshold
        self._groups = []
        self._max_dens = {}
        mylog.info("Initializing HOP")
        self.__obtain_particles()
        self.__run_hop()
        mylog.info("Parsing outputs")
        self.__parse_output()
        mylog.debug("Finished. (%s)", len(self))

    def __obtain_particles(self):
        if self.dm_only: ii = self.__get_dm_indices()
        else: ii = slice(None)
        self.particle_fields = {}
        for field in ["particle_position_%s" % ax for ax in 'xyz'] + \
                     ["ParticleMassMsun"]:
            tot_part = self.data_source[field].size
            self.particle_fields[field] = self.data_source[field][ii]
        self._base_indices = na.arange(tot_part)[ii]

    def __run_hop(self):
        self.densities, self.tags = \
            RunHOP(self.particle_fields["particle_position_x"],
                   self.particle_fields["particle_position_y"],
                   self.particle_fields["particle_position_z"],
                   self.particle_fields["ParticleMassMsun"],
                   self.threshold)
        self.particle_fields["densities"] = self.densities
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
        dens = self.densities[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i+1]
            if i == -1:
                cp += counts[i+1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(HopGroup(self, i, group_indices))
            md_i = na.argmax(dens[cp:cp_c])
            px, py, pz = [self.particle_fields['particle_position_%s'%ax][group_indices]
                                            for ax in 'xyz']
            self._max_dens[i] = (dens[cp:cp_c][md_i],
                                 px[md_i], py[md_i], pz[md_i])
            cp += counts[i+1]

    def __len__(self):
        return len(self._groups)
 
    def __iter__(self):
        return HopIterator(self)

    def __getitem__(self, key):
        return self._groups[key]

    def write_out(self, filename="HopAnalysis.out"):
        """
        Write out standard HOP information to *filename*.
        """
        f = open(filename,"w")
        f.write("\t".join(["# Group","Mass","# part","max dens"
                           "x","y","z", "center-of-mass",
                           "x","y","z",
                           "vx","vy","vz","max_r","\n"]))
        for group in self:
            f.write("%10i\t" % group.id)
            f.write("%0.9e\t" % group.total_mass())
            f.write("%10i\t" % group.indices.size)
            f.write("%0.9e\t" % group.maximum_density())
            f.write("\t".join(["%0.9e" % v for v in group.maximum_density_location()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.center_of_mass()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.bulk_velocity()]))
            f.write("\t")
            f.write("%0.9e\t" % group.maximum_radius())
            f.write("\n")
        f.close()

class HopIterator(object):
    def __init__(self, hop):
        self.hop = hop
        self.index = -1

    def next(self):
        self.index += 1
        if self.index == len(self.hop): raise StopIteration
        return self.hop[self.index]

class HopGroup(object):
    """
    A data source that returns particle information about the members of a
    HOP-identified halo.
    """
    def __init__(self, hop_output, id, indices):
        self.hop_output = hop_output
        self.id = id
        self.data = hop_output.data_source
        self.indices = hop_output._base_indices[indices]
        
    def center_of_mass(self):
        """
        Calculate and return the center of mass.
        """
        c_vec = self.maximum_density_location() - na.array([0.5,0.5,0.5])
        pm = self["ParticleMassMsun"]
        cx = (self["particle_position_x"] - c_vec[0])
        cy = (self["particle_position_y"] - c_vec[1])
        cz = (self["particle_position_z"] - c_vec[2])
        com = na.array([v-na.floor(v) for v in [cx,cy,cz]])
        return (com*pm).sum(axis=1)/pm.sum() + c_vec

    def maximum_density(self):
        """
        Return the HOP-identified maximum density.
        """
        return self.hop_output._max_dens[self.id][0]

    def maximum_density_location(self):
        """
        Return the location HOP identified as maximally dense.
        """
        return na.array([
                self.hop_output._max_dens[self.id][1],
                self.hop_output._max_dens[self.id][2],
                self.hop_output._max_dens[self.id][3]])

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

    def maximum_radius(self, center_of_mass=True):
        """
        Returns the maximum radius in the halo for all particles,
        either from the point of maximum density or from the (default)
        *center_of_mass*.
        """
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        rx = na.abs(self["particle_position_x"]-center[0])
        ry = na.abs(self["particle_position_y"]-center[1])
        rz = na.abs(self["particle_position_z"]-center[2])
        r = na.sqrt(na.minimum(rx, 1.0-rx)**2.0
                +   na.minimum(ry, 1.0-ry)**2.0
                +   na.minimum(rz, 1.0-rz)**2.0)
        return r.max()

    def __getitem__(self, key):
        return self.data[key][self.indices]

    def get_sphere(self, center_of_mass=True):
        """
        Returns an EnzoSphere centered on either the point of maximum density
        or the *center_of_mass*, with the maximum radius of the halo.
        """
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        radius = self.maximum_radius()
        # A bit of a long-reach here...
        sphere = self.hop_output.data_source.hierarchy.sphere(
                        center, radius=radius)
        return sphere

class HaloFinder(ParallelAnalysisInterface):
    def __init__(self, pf, threshold=160.0, dm_only=True):
        self.pf = pf
        self.hierarchy = pf.hierarchy
        # do it once with no padding so the total_mass is correct (no duplicated particles)
        self.padding = 0.0
        LE, RE, self.source = self._partition_hierarchy_3d(padding=self.padding)
        # For scaling the threshold, note that it's a passthrough
        total_mass = self._mpi_allsum(self.source["ParticleMassMsun"].sum())
        self.padding = 0.2 #* pf["unitary"]
        LE, RE, self.source = self._partition_hierarchy_3d(padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        self._reposition_particles((LE, RE))
        # calculate hop on each sub region
        hop_list = HopList(self.source, threshold, dm_only)
        # include only haloes that reside in the real part of the box
        self._parse_hoplist(hop_list)
        # collect all the haloes into one big ordered list
        self._join_hoplists(hop_list)

    @parallel_passthrough
    # on each processor, go through the hoplist and throw away haloes that aren't in the 'real'
    # part of the subbox, and add the good ones to a list, _groups
    def _parse_hoplist(self, hop_list):
        LE, RE = bounds
        # find the unique halo ids
        unique_ids = na.unique(hop_list.tags)
        # find the size of each halo
        counts = na.bincount(hop_list.tags)
        # sort the particles by their halo ID
        sort_indices = na.argsort(hop_list.tags)
        # grab the particles into a simple list
        grab_indices = na.indices(hop_list.tags.shape).ravel()
        # get the densities of the particles, ordered as above
        dens = hop_list.densities[sort_indices]
        cp = 0
        print 'uids %d' % (unique_ids.size)
        # now we'll loop over the haloes, and only keep the real ones
        for ii in unique_ids:
            cp_c = cp + counts[ii+1]
            if ii == -1:
                cp += counts[ii+1]
                continue
            group_indices = grab_indices[cp:cp_c]
            # find the id and position of the densest particle
            md_i = na.argmax(dens[cp:cp_c])
            px = self.particle_fields["particle_position_x"][group_indices]
            px, py, pz = [self.particle_fields["particle_position_%s"%ax][group_indices] for ax in 'xyz']
            max_dens = (dens[cp:cp_c][md_i], px[md_i], py[md_i], pz[md_i])
            # if the most dense particle is in the box, keep it
            if ((max_dens[1:3] >= LE+self.padding) && (max_dens[1:3] < RE-self.padding)):
                self._groups.append(HopGroup(self, ii, group_indices)
                self._max_dens[ii] = max_dens
            cp += counts[ii+1]
    
    @parallel_passthrough
    def _join_hoplists(self, hop_list):
        # First we get the total number of halos the entire collection
        # has identified
        nhalos = self._mpi_allsum(len(_groups))
        # I don't know how this stuff works, so I'm going to sketch it out.
        # collect all the HopGroups into one big list
        haloes = self._mpi_allcollectwhatever(_groups)
        # sort the list by the size of the groups
        haloes.sort(lambda x, y: cmp(len(x.indices),len(y.indices)))
        # reassign their ID
        for i,halo in enumerate(haloes):
            halo.id = i
        
    
    @parallel_passthrough
    def _reposition_particles(self, bounds):
        # This only does periodicity.  We do NOT want to deal with anything
        # else.  The only reason we even do periodicity is the 
        LE, RE = bounds
        dw = self.pf["DomainRightEdge"] - self.pf["DomainLeftEdge"]
        for i, ax in enumerate('xyz'):
            arr = self.source["particle_position_%s" % ax]
            arr[arr < LE[i]-self.padding] += dw[i]
            arr[arr > RE[i]+self.padding] -= dw[i]

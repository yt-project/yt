"""
HOP-output data handling

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
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
        self.data_source = data_source
        self.dm_only = dm_only
        self.threshold = threshold
        self._groups = []
        self._max_dens = {}
        mylog.info("Initializing HOP")
        self.__run_hop()
        mylog.info("Parsing outputs")
        self.__parse_output()
        mylog.debug("Finished.")

    def __run_hop(self):
        if self.dm_only: ii = self.__get_dm_indices()
        else: ii = slice(None)
        self._base_indices = ii
        self.densities, self.tags = \
            RunHOP(self.data_source["particle_position_x"][ii],
                   self.data_source["particle_position_y"][ii],
                   self.data_source["particle_position_z"][ii],
                   self.data_source["ParticleMassMsun"][ii],
                   self.data_source["particle_index"][ii].astype('int64'),
                   self.threshold)

    def __get_dm_indices(self):
        if 'creation_time' in self.data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on creation time")
            return (self.data_source["creation_time"] < 0)
        elif 'particle_type' in self.data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on particle type")
            return (self.data_source["particle_type"] == 1)
        else:
            raise KeyError

    def __parse_output(self):
        unique_ids = na.unique(self.tags)
        counts = na.bincount((self.tags+1)
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i+1]
            if i == -1:
                cp += counts[i+1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(HopGroup(self, i, group_indices))
            md_i = na.argmax(self.densities[sort_indices][cp:cp_c])
            px, py, pz = [self.data_source['particle_position_%s'%ax][self._base_indices][group_indices]
                                            for ax in 'xyz']
            self._max_dens[i] = (self.densities[sort_indices][cp:cp_c][md_i],
                                 px[md_i], py[md_i], pz[md_i])
            cp += counts[i+1]

    def __len__(self):
        return len(self._groups)
 
    def __iter__(self):
        return HopIterator(self)

    def __getitem__(self, key):
        return self._groups[key]

    def write_out(self, filename="HopAnalysis.out"):
        f = open(filename,"w")
        f.write("# Center of mass does NOT account for periodicity!\n")
        f.write("\t".join(["# Group","Mass","# part","max dens"
                           "x","y","z", "center-of-mass",
                           "x","y","z",
                           "vx","vy","vz"]))
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

    def __init__(self, hop_output, id, indices):
        self.hop_output = hop_output
        self.id = id
        self.data = hop_output.data_source
        self.indices = indices
        self._base_indices = hop_output._base_indices
        
    def center_of_mass(self):
        # Center of mass does not account for periodicity and is likely wrong!
        pm = self["ParticleMassMsun"]
        cx = (self["particle_position_x"] * pm).sum()
        cy = (self["particle_position_y"] * pm).sum()
        cz = (self["particle_position_z"] * pm).sum()
        return na.array([cx,cy,cz])/pm.sum()

    def maximum_density(self):
        return self.hop_output._max_dens[self.id][0]

    def maximum_density_location(self):
        return (self.hop_output._max_dens[self.id][1],
                self.hop_output._max_dens[self.id][2],
                self.hop_output._max_dens[self.id][3])

    def total_mass(self):
        return self["ParticleMassMsun"].sum()

    def bulk_velocity(self):
        pm = self["ParticleMassMsun"]
        vx = (self["particle_velocity_x"] * pm).sum()
        vy = (self["particle_velocity_y"] * pm).sum()
        vz = (self["particle_velocity_z"] * pm).sum()
        return na.array([vx,vy,vz])/pm.sum()

    def maximum_radius(self):
        center = self.maximum_density_location()
        #center = self.center_of_mass()
        rx = na.abs(self["particle_position_x"]-center[0])
        ry = na.abs(self["particle_position_y"]-center[1])
        rz = na.abs(self["particle_position_z"]-center[2])
        r = na.sqrt(na.minimum(rx, 1.0-rx)**2.0
                +   na.minimum(ry, 1.0-ry)**2.0
                +   na.minimum(rz, 1.0-rz)**2.0)
        return r.max()

    def __getitem__(self, key):
        return self.data[key][self._base_indices][self.indices]

    def get_sphere(self, center_of_mass=True):
        if center_of_mass: center = self.center_of_mass()
        else: center = self.maximum_density_location()
        radius = self.maximum_radius()
        # A bit of a long-reach here...
        sphere = self.hop_output.data_source.hierarchy.sphere(
                        center, radius=radius)
        return sphere

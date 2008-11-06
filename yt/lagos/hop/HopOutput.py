"""
HOP-output data handling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
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
from numpy import *

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
        self._groups2 = {}
        self._max_dens = {}
        self._max_dens2 = {}
        mylog.info("Initializing HOP")
        self.__obtain_particles()
        self.__enlarge_data()
        self.__slice_data()
        self.__run_hops()
        self.__reduce_data()
        #self.__run_hop()
        mylog.info("Parsing outputs")
        self.__parse_outputs()
        #self.__parse_output()
        self.__glue_outputs()
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

    def __enlarge_data(self):
        sizetemp = self.particle_fields["particle_position_x"].size
        self.tempx = [i for i in range(3*sizetemp)]
        self.tempy = [i for i in range(3*sizetemp)]
        self.tempz = [i for i in range(3*sizetemp)]
        self.tempm = [i for i in range(3*sizetemp)]
        self.tracking = [i for i in range(3*sizetemp)]
        # reverse below 0, straight up between 0 and 1, reverse above 1
        for part,i in enumerate(self.particle_fields["particle_position_x"]):
            self.tempx[sizetemp - 1 - part] = -i
            self.tempx[sizetemp + part] = i
            self.tempx[3*sizetemp - 1 - part] = 2 - i
            self.tracking[sizetemp - 1 - part] = part
            self.tracking[sizetemp + part] = part
            self.tracking[3*sizetemp - 1 - part] = part
        for part,i in enumerate(self.particle_fields["particle_position_y"]):
            self.tempy[sizetemp - 1 - part] = -i
            self.tempy[sizetemp + part] = i
            self.tempy[3*sizetemp - 1 - part] = 2 - i
        for part,i in enumerate(self.particle_fields["particle_position_z"]):
            self.tempz[sizetemp - 1 - part] = -i
            self.tempz[sizetemp + part] = i
            self.tempz[3*sizetemp - 1 - part] = 2 - i
        for part,i in enumerate(self.particle_fields["ParticleMassMsun"]):
            self.tempm[sizetemp - 1 - part] = i
            self.tempm[sizetemp + part] = i
            self.tempm[3*sizetemp -1 - part] = i

    def __slice_data(self):
        cuts = 3
        padding = 0.2
        self.temp2x = {}
        self.temp2y = {}
        self.temp2z = {}
        self.temp2m = {}
        self.tracking2 = {}
        # loop over the sub-boxes
        for i in range(cuts):
            for j in range(cuts):
                for k in range(cuts):
                    x_t = []
                    y_t = []
                    z_t = []
                    m_t = []
                    t_t = {}
                    # define the edges
                    xleft = (1.0 * i)/cuts - padding
                    xright = (1.0 * (i+1))/cuts + padding
                    yleft = (1.0 * j)/cuts - padding
                    yright = (1.0 * (j+1))/cuts + padding
                    zleft = (1.0 * k)/cuts - padding
                    zright = (1.0 * (k+1))/cuts + padding
                    # loop over the temporary expanded fields, test to see if they're in the sub-box
                    # and then re-map to a unit sized box, we'll have inclusive boundary on the left,
                    # exclusive on the right
                    jj = 0
                    for part,ii in enumerate(self.tempx):
                        if ((ii >= xleft) and (ii < xright) and (self.tempy[part] >= yleft) and \
                        (self.tempy[part] < yright) and (self.tempz[part] >= zleft) and \
                        (self.tempz[part] < zright)):
                            x_t.append((ii-xleft)*(1.0/(xright-xleft)))
                            y_t.append((self.tempy[part]-yleft)*(1.0/(yright-yleft)))
                            z_t.append((self.tempz[part]-zleft)*(1.0/(zright-zleft)))
                            m_t.append(self.tempm[part])
                            t_t[jj] = self.tracking[part] # put real ID in t_t at jj
                            jj += 1
                    # save them to their position in the dict
                    self.temp2x[((i*cuts) + j)*cuts + k] = array(x_t)
                    self.temp2y[((i*cuts) + j)*cuts + k] = array(y_t)
                    self.temp2z[((i*cuts) + j)*cuts + k] = array(z_t)
                    self.temp2m[((i*cuts) + j)*cuts + k] = array(m_t)
                    self.tracking2[((i*cuts) + j)*cuts + k] = t_t

    def __reduce_data(self):
        cuts = 3
        padding = 0.2
        # loop over the sub-boxes
        for i in range(cuts):
            for j in range(cuts):
                for k in range(cuts):
                    xleft = (1.0 * i)/cuts - padding
                    xright = (1.0 * (i+1))/cuts + padding
                    yleft = (1.0 * j)/cuts - padding
                    yright = (1.0 * (j+1))/cuts + padding
                    zleft = (1.0 * k)/cuts - padding
                    zright = (1.0 * (k+1))/cuts + padding
                    # we're going to reduce the values of the particles to their true position
                    # with some going negative, of course
                    self.temp2x[((i*cuts) + j)*cuts + k] = self.temp2x[((i*cuts) + j)*cuts + k]*(xright-xleft)+xleft
                    self.temp2y[((i*cuts) + j)*cuts + k] = self.temp2y[((i*cuts) + j)*cuts + k]*(yright-yleft)+yleft
                    self.temp2z[((i*cuts) + j)*cuts + k] = self.temp2z[((i*cuts) + j)*cuts + k]*(zright-zleft)+zleft

    def __run_hops(self):
        cuts = 3
        padding = 0.2
        nParts = self.particle_fields["particle_position_x"].size
        dens_fix = 1. / (2*padding + 1./cuts)**3
        self.densities2 = {}
        self.tags2 = {}
        for i in range(cuts):
            for j in range(cuts):
                for k in range(cuts):
                    print "run_hops %d" % (((i*cuts) + j)*cuts + k)
                    size = float(self.temp2x[((i*cuts) + j)*cuts + k].size)
                    print "size %d" % size
                    self.densities2[((i*cuts) + j)*cuts + k], self.tags2[((i*cuts) + j)*cuts + k] = \
                        RunHOP(self.temp2x[((i*cuts) + j)*cuts + k],
                               self.temp2y[((i*cuts) + j)*cuts + k],
                               self.temp2z[((i*cuts) + j)*cuts + k],
                               self.temp2m[((i*cuts) + j)*cuts + k],
                               self.threshold,size/nParts*dens_fix)

    def __run_hop(self):
        self.densities, self.tags = \
            RunHOP(self.particle_fields["particle_position_x"],
                   self.particle_fields["particle_position_y"],
                   self.particle_fields["particle_position_z"],
                   self.particle_fields["ParticleMassMsun"],
                   self.threshold, 1.0)
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

    def __parse_outputs(self):
        cuts = 3
        for i in range(cuts):
            for j in range(cuts):
                for k in range(cuts):
                    self._groups2[((i*cuts) + j)*cuts + k] = []
                    self._max_dens2[((i*cuts) + j)*cuts + k] = []
                    unique_ids = na.unique(self.tags2[((i*cuts) + j)*cuts + k])
                    counts = na.bincount(self.tags2[((i*cuts) + j)*cuts + k]+1)
                    sort_indices = na.argsort(self.tags2[((i*cuts) + j)*cuts + k])
                    grab_indices = na.indices(self.tags2[((i*cuts) + j)*cuts + k].shape).ravel()[sort_indices]
                    dens = self.densities2[((i*cuts) + j)*cuts + k][sort_indices]
                    cp = 0
                    print 'uids %d' % (unique_ids.size)
                    for ii in unique_ids:
                        cp_c = cp + counts[ii+1]
                        if ii == -1:
                            cp += counts[ii+1]
                            continue
                        group_indices = grab_indices[cp:cp_c]
                        md_i = na.argmax(dens[cp:cp_c])
                        px = self.temp2x[((i*cuts) + j)*cuts + k][group_indices]
                        py = self.temp2y[((i*cuts) + j)*cuts + k][group_indices]
                        pz = self.temp2z[((i*cuts) + j)*cuts + k][group_indices]
                        max_dens = (dens[cp:cp_c][md_i], px[md_i], py[md_i], pz[md_i])
                        # as a shortcut, if the most dense particle is in the current box, we keep the group
                        if ((math.floor(max_dens[1]*cuts) == i) and (math.floor(max_dens[2]*cuts) == j) \
                                and (math.floor(max_dens[3]*cuts) == k)):
                            # I think this will make the group_indices correct in the overall 'real' sense
                            temp_real_indices = arange(group_indices.size)
                            for gg,gi in enumerate(group_indices):
                                temp_real_indices[gg] = self.tracking2[((i*cuts) + j)*cuts + k][gi]
                            self._groups2[((i*cuts) + j)*cuts + k].append(HopGroup(self,ii, temp_real_indices))
                            self._max_dens2[((i*cuts) + j)*cuts + k].append(max_dens)
                        cp += counts[ii+1]
        for i in range(cuts**3):
            print '%d groups2 %d' % (i,len(self._groups2[i]))

    def __glue_outputs(self):
        cuts = 3
        ii = 0
        iii = 0
        for i in range(cuts):
            for j in range(cuts):
                for k in range(cuts):
                    print '%d len = %d' % (((i*cuts) + j)*cuts + k,len(self._groups2[((i*cuts) + j)*cuts + k]))
                    for group2 in self._groups2[((i*cuts) + j)*cuts + k]:
                        self._groups.append(HopGroup(self,iii,group2.indices))
                        iii += 1
                    for dens in self._max_dens2[((i*cuts) + j)*cuts + k]:
                        self._max_dens[ii] = dens
                        ii += 1

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

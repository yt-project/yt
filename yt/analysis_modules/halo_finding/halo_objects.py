"""
HOP-output data handling

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Stephen Skory <s@skory.us>
Affiliation: UCSD Physics/CASS
Author: Geoffrey So <gsiisg@gmail.com> (Ellipsoidal functions)
Affiliation: UCSD Physics/CASS
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

import gc
import h5py
import itertools
import math
import numpy as na
import random
import sys
import os.path as path
from collections import defaultdict

from yt.funcs import *

from yt.config import ytcfg
from yt.utilities.performance_counters import \
    yt_counters, time_function
from yt.utilities.math_utils import periodic_dist
from yt.utilities.physical_constants import rho_crit_now, mass_sun_cgs

from .hop.EnzoHop import RunHOP
from .fof.EnzoFOF import RunFOF
try:
    from parallel_hop.parallel_hop_interface import \
        ParallelHOPHaloFinder
except ImportError:
    mylog.debug("Parallel HOP not imported.")

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelDummy, \
    ParallelAnalysisInterface, \
    parallel_blocking_call

TINY = 1.e-40

# Ellipsoid funtions.
# define the rotation matrix needed later
def RX(ax):
    rot_matrix = na.array([[1, 0, 0], [0, na.cos(ax), na.sin(ax)],
        [0, -na.sin(ax), na.cos(ax)]])
    return rot_matrix
def RY(ay):
    rot_matrix = na.array([[na.cos(ay), 0, -na.sin(ay)], [0, 1, 0],
        [na.sin(ay), 0, na.cos(ay)]])
    return rot_matrix
def RZ(az):
    rot_matrix = na.array([[na.cos(az), na.sin(az), 0],
        [-na.sin(az), na.cos(az), 0], [0, 0, 1]])
    return rot_matrix

class Halo(object):
    """
    A data source that returns particle information about the members of a
    HOP-identified halo.
    """
    __metaclass__ = ParallelDummy  # This will proxy up our methods
    _distributed = False
    _processing = False
    _owner = 0
    indices = None
    dont_wrap = ["get_sphere", "write_particle_list"]
    extra_wrap = ["__getitem__"]

    def __init__(self, halo_list, id, indices=None, size=None, CoM=None,
        max_dens_point=None, group_total_mass=None, max_radius=None,
        bulk_vel=None, tasks=None, rms_vel=None):
        self._max_dens = halo_list._max_dens
        self.id = id
        self.data = halo_list._data_source
        self.pf = self.data.pf
        self.gridsize = (self.pf.domain_right_edge - \
                 self.pf.domain_left_edge)
        if indices is not None:
            self.indices = halo_list._base_indices[indices]
        else:
            self.indices = None
        # We assume that if indices = None, the instantiator has OTHER plans
        # for us -- i.e., setting it somehow else
        self.size = size
        self.CoM = CoM
        self.max_dens_point = max_dens_point
        self.group_total_mass = group_total_mass
        self.max_radius = max_radius
        self.bulk_vel = bulk_vel
        self.tasks = tasks
        self.rms_vel = rms_vel
        self.bin_count = None
        self.overdensity = None

    def center_of_mass(self):
        r"""Calculate and return the center of mass.

        The center of mass of the halo is directly calculated and returned.

        Examples
        --------
        >>> com = halos[0].center_of_mass()
        """
        if self.CoM is not None:
            return self.CoM
        pm = self["ParticleMassMsun"]
        cx = self["particle_position_x"]
        cy = self["particle_position_y"]
        cz = self["particle_position_z"]
        if isinstance(self, FOFHalo):
            c_vec = na.array([cx[0], cy[0], cz[0]]) - self.pf.domain_center
        else:
            c_vec = self.maximum_density_location() - self.pf.domain_center
        cx = (cx - c_vec[0])
        cy = (cy - c_vec[1])
        cz = (cz - c_vec[2])
        com = na.array([v - na.floor(v) for v in [cx, cy, cz]])
        return (com * pm).sum(axis=1) / pm.sum() + c_vec

    def maximum_density(self):
        r"""Return the HOP-identified maximum density. Not applicable to
        FOF halos.

        Return the HOP-identified maximum density. Not applicable to FOF halos.

        Examples
        --------
        >>> max_dens = halos[0].maximum_density()
        """
        if self.max_dens_point is not None:
            return self.max_dens_point[0]
        return self._max_dens[self.id][0]

    def maximum_density_location(self):
        r"""Return the location HOP identified as maximally dense. Not
        applicable to FOF halos.

        Return the location HOP identified as maximally dense.

        Examples
        --------
        >>> max_dens_loc = halos[0].maximum_density_location()
        """
        if self.max_dens_point is not None:
            return self.max_dens_point[1:]
        return na.array([
                self._max_dens[self.id][1],
                self._max_dens[self.id][2],
                self._max_dens[self.id][3]])

    def total_mass(self):
        r"""Returns the total mass in solar masses of the halo.

        Returns the total mass in solar masses of just the particles in the
        halo.

        Examples
        --------
        >>> halos[0].total_mass()
        """
        if self.group_total_mass is not None:
            return self.group_total_mass
        return self["ParticleMassMsun"].sum()

    def bulk_velocity(self):
        r"""Returns the mass-weighted average velocity in cm/s.

        This calculates and returns the mass-weighted average velocity of just
        the particles in the halo in cm/s.

        Examples
        --------
        >>> bv = halos[0].bulk_velocity()
        """
        if self.bulk_vel is not None:
            return self.bulk_vel
        pm = self["ParticleMassMsun"]
        vx = (self["particle_velocity_x"] * pm).sum()
        vy = (self["particle_velocity_y"] * pm).sum()
        vz = (self["particle_velocity_z"] * pm).sum()
        return na.array([vx, vy, vz]) / pm.sum()

    def rms_velocity(self):
        r"""Returns the mass-weighted RMS velocity for the halo
        particles in cgs units.

        Calculate and return the mass-weighted RMS velocity for just the
        particles in the halo.  The bulk velocity of the halo is subtracted
        before computation.

        Examples
        --------
        >>> rms_vel = halos[0].rms_velocity()
        """
        if self.rms_vel is not None:
            return self.rms_vel
        bv = self.bulk_velocity()
        pm = self["ParticleMassMsun"]
        sm = pm.sum()
        vx = (self["particle_velocity_x"] - bv[0]) * pm / sm
        vy = (self["particle_velocity_y"] - bv[1]) * pm / sm
        vz = (self["particle_velocity_z"] - bv[2]) * pm / sm
        s = vx ** 2. + vy ** 2. + vz ** 2.
        ms = na.mean(s)
        return na.sqrt(ms) * pm.size

    def maximum_radius(self, center_of_mass=True):
        r"""Returns the maximum radius in the halo for all particles,
        either from the point of maximum density or from the
        center of mass.

        The maximum radius from the most dense point is calculated.  This
        accounts for periodicity.

        Parameters
        ----------
        center_of_mass : bool
            True chooses the center of mass when
            calculating the maximum radius.
            False chooses from the maximum density location for HOP halos
            (it has no effect for FOF halos).
            Default = True.

        Examples
        --------
        >>> radius = halos[0].maximum_radius()
        """
        if self.max_radius is not None:
            return self.max_radius
        if center_of_mass:
            center = self.center_of_mass()
        else:
            center = self.maximum_density_location()
        rx = na.abs(self["particle_position_x"] - center[0])
        ry = na.abs(self["particle_position_y"] - center[1])
        rz = na.abs(self["particle_position_z"] - center[2])
        DW = self.data.pf.domain_right_edge - self.data.pf.domain_left_edge
        r = na.sqrt(na.minimum(rx, DW[0] - rx) ** 2.0
                + na.minimum(ry, DW[1] - ry) ** 2.0
                + na.minimum(rz, DW[2] - rz) ** 2.0)
        return r.max()

    def __getitem__(self, key):
        if ytcfg.getboolean("yt", "inline") == False:
            return self.data[key][self.indices]
        else:
            return self.data[key][self.indices]

    def get_sphere(self, center_of_mass=True):
        r"""Returns a sphere source.

        This will generate a new, empty sphere source centered on this halo,
        with the maximum radius of the halo. This can be used like any other
        data container in yt.

        Parameters
        ----------
        center_of_mass : bool, optional
            True chooses the center of mass when
            calculating the maximum radius.
            False chooses from the maximum density location for HOP halos
            (it has no effect for FOF halos).
            Default = True.

        Returns
        -------
        sphere : `yt.data_objects.api.AMRSphereBase`
            The empty data source.

        Examples
        --------
        >>> sp = halos[0].get_sphere()
        """
        if center_of_mass:
            center = self.center_of_mass()
        else:
            center = self.maximum_density_location()
        radius = self.maximum_radius()
        # A bit of a long-reach here...
        sphere = self.data.hierarchy.sphere(
                        center, radius=radius)
        return sphere

    def get_size(self):
        if self.size is not None:
            return self.size
        return self.indices.size

    def write_particle_list(self, handle):
        self._processing = True
        gn = "Halo%08i" % (self.id)
        handle.create_group("/%s" % gn)
        for field in ["particle_position_%s" % ax for ax in 'xyz'] \
                   + ["particle_velocity_%s" % ax for ax in 'xyz'] \
                   + ["particle_index"] + ["ParticleMassMsun"]:
            handle.create_dataset("/%s/%s" % (gn, field), data=self[field])
        if 'creation_time' in self.data.pf.h.field_list:
            handle.create_dataset("/%s/creation_time" % gn,
                data=self['creation_time'])
        n = handle["/%s" % gn]
        # set attributes on n
        self._processing = False

    def virial_mass(self, virial_overdensity=200., bins=300):
        r"""Return the virial mass of the halo in Msun,
        using only the particles
        in the halo (no baryonic information used).

        The virial mass is calculated, using the built in `Halo.virial_info`
        functionality.  The mass is then returned.

        Parameters
        ----------
        virial_overdensity : float
            The overdensity threshold compared to the universal average when
            calculating the virial mass. Default = 200.
        bins : int
            The number of spherical bins used to calculate overdensities.
            Default = 300.

        Returns
        -------
        mass : float
            The virial mass in solar masses of the particles in the halo.  -1
            if not virialized.

        Examples
        --------
        >>> vm = halos[0].virial_mass()
        """
        self.virial_info(bins=bins)
        vir_bin = self.virial_bin(virial_overdensity=virial_overdensity,
            bins=bins)
        if vir_bin != -1:
            return self.mass_bins[vir_bin]
        else:
            return -1

    def virial_radius(self, virial_overdensity=200., bins=300):
        r"""Return the virial radius of the halo in code units.

        The virial radius of the halo is calculated, using only the particles
        in the halo (no baryonic information used). Returns -1 if the halo is
        not virialized.

        Parameters
        ----------
        virial_overdensity : float
            The overdensity threshold compared to the universal average when
            calculating the virial radius. Default = 200.
        bins : integer
            The number of spherical bins used to calculate overdensities.
            Default = 300.

        Returns
        -------
        radius : float
            The virial raius in code units of the particles in the halo.  -1
            if not virialized.

        Examples
        --------
        >>> vr = halos[0].virial_radius()
        """
        self.virial_info(bins=bins)
        vir_bin = self.virial_bin(virial_overdensity=virial_overdensity,
            bins=bins)
        if vir_bin != -1:
            return self.radial_bins[vir_bin]
        else:
            return -1

    def virial_bin(self, virial_overdensity=200., bins=300):
        r"""Returns the bin index of the virial radius of the halo. Generally,
        it is better to call virial_radius instead, which calls this function
        automatically.
        """
        self.virial_info(bins=bins)
        over = (self.overdensity > virial_overdensity)
        if (over == True).any():
            vir_bin = max(na.arange(bins + 1)[over])
            return vir_bin
        else:
            return -1

    def virial_info(self, bins=300):
        r"""Calculates the virial information for the halo. Generally, it is
        better to call virial_radius or virial_mass instead, which calls this
        function automatically.
        """
        # Skip if we've already calculated for this number of bins.
        if self.bin_count == bins and self.overdensity is not None:
            return None
        self.bin_count = bins
        # Cosmology
        h = self.pf.hubble_constant
        Om_matter = self.pf.omega_matter
        z = self.pf.current_redshift
        period = self.pf.domain_right_edge - \
            self.pf.domain_left_edge
        cm = self.pf["cm"]
        thissize = max(self.size, self.indices.size)
        rho_crit = rho_crit_now * h ** 2.0 * Om_matter  # g cm^-3
        Msun2g = mass_sun_cgs
        rho_crit = rho_crit * ((1.0 + z) ** 3.0)
        # Get some pertinent information about the halo.
        self.mass_bins = na.zeros(self.bin_count + 1, dtype='float64')
        dist = na.empty(thissize, dtype='float64')
        cen = self.center_of_mass()
        mark = 0
        # Find the distances to the particles. I don't like this much, but I
        # can't see a way to eliminate a loop like this, either here or in
        # yt.math.
        for pos in itertools.izip(self["particle_position_x"],
                self["particle_position_y"], self["particle_position_z"]):
            dist[mark] = periodic_dist(cen, pos, period)
            mark += 1
        # Set up the radial bins.
        # Multiply min and max to prevent issues with digitize below.
        self.radial_bins = na.logspace(math.log10(min(dist) * .99 + TINY),
            math.log10(max(dist) * 1.01 + 2 * TINY), num=self.bin_count + 1)
        # Find out which bin each particle goes into, and add the particle
        # mass to that bin.
        inds = na.digitize(dist, self.radial_bins) - 1
        if self["particle_position_x"].size > 1:
            for index in na.unique(inds):
                self.mass_bins[index] += \
                na.sum(self["ParticleMassMsun"][inds == index])
        # Now forward sum the masses in the bins.
        for i in xrange(self.bin_count):
            self.mass_bins[i + 1] += self.mass_bins[i]
        # Calculate the over densities in the bins.
        self.overdensity = self.mass_bins * Msun2g / \
        (4./3. * math.pi * rho_crit * \
        (self.radial_bins * cm)**3.0)
        
    def _get_ellipsoid_parameters_basic(self):
        na.seterr(all='ignore')
        # check if there are 4 particles to form an ellipsoid
        # neglecting to check if the 4 particles in the same plane,
        # that is almost certainly never to occur,
        # will deal with it later if it ever comes up
        if na.size(self["particle_position_x"]) < 4:
            print "not enough particles to form ellipsoid returning zeros"
            return (0, 0, 0, 0, 0, 0, 0)
        # Calculate the parameters that describe the ellipsoid of
        # the particles that constitute the halo. This function returns
        # all the parameters except for the center of mass.
        com = self.center_of_mass()
        position = [self["particle_position_x"],
		    self["particle_position_y"],
		    self["particle_position_z"]]
        # Locate the furthest particle from com, its vector length and index
	DW = na.array([self.gridsize[0],self.gridsize[1],self.gridsize[2]])
	position = [position[0] - com[0],
		    position[1] - com[1],
		    position[2] - com[2]]
	# different cases of particles being on other side of boundary
	for axis in range(na.size(DW)):
	    cases = na.array([position[axis],
	  		      position[axis] + DW[axis],
			      position[axis] - DW[axis]])        
            # pick out the smallest absolute distance from com
            position[axis] = na.choose(na.abs(cases).argmin(axis=0), cases)
	# find the furthest particle's index
	r = na.sqrt(position[0]**2 +
		    position[1]**2 +
		    position[2]**2)
        A_index = r.argmax()
        mag_A = r.max()
        # designate the A vector
	A_vector = (position[0][A_index],
		    position[1][A_index],
		    position[2][A_index])
        # designate the e0 unit vector
        e0_vector = A_vector / mag_A
        # locate the tB particle position by finding the max B
	e0_vector_copy = na.empty((na.size(position[0]), 3), dtype='float64')
        for i in range(3):
            e0_vector_copy[:, i] = e0_vector[i]
        rr = na.array([position[0],
		       position[1],
		       position[2]]).T # Similar to tB_vector in old code.
        tC_vector = na.cross(e0_vector_copy, rr)
        te2 = tC_vector.copy()
        for dim in range(3):
            te2[:,dim] *= na.sum(tC_vector**2., axis = 1)**(-0.5)
        te1 = na.cross(te2, e0_vector_copy)
        length = na.abs(-na.sum(rr * te1, axis = 1) * \
            (1. - na.sum(rr * e0_vector_copy, axis = 1)**2. * \
            mag_A**-2.)**(-0.5))
        # This problem apparently happens sometimes, that the NaNs are turned
        # into infs, which messes up the nanargmax below.
        length[length == na.inf] = 0.
        tB_index = na.nanargmax(length) # ignores NaNs created above.
        mag_B = length[tB_index]
        e1_vector = te1[tB_index]
        e2_vector = te2[tB_index]
        temp_e0 = rr.copy()
        temp_e1 = rr.copy()
        temp_e2 = rr.copy()
        for dim in range(3):
            temp_e0[:,dim] = e0_vector[dim]
            temp_e1[:,dim] = e1_vector[dim]
            temp_e2[:,dim] = e2_vector[dim]
        length = na.abs(na.sum(rr * temp_e2, axis = 1) * (1 - \
            na.sum(rr * temp_e0, axis = 1)**2. * mag_A**-2. - \
            na.sum(rr * temp_e1, axis = 1)**2. * mag_B**-2)**(-0.5))
        length[length == na.inf] = 0.
        tC_index = na.nanargmax(length)
        mag_C = length[tC_index]
        # tilt is calculated from the rotation about x axis
        # needed to align e1 vector with the y axis
        # after e0 is aligned with x axis
        # find the t1 angle needed to rotate about z axis to align e0 to x
        t1 = na.arctan(e0_vector[1] / e0_vector[0])
        r1 = (e0_vector * RZ(t1).transpose()).sum(axis = 1)
        # find the t2 angle needed to rotate about y axis to align e0 to x
        t2 = na.arctan(-r1[2] / r1[0])
        r2 = na.dot(RY(t2), na.dot(RZ(t1), e1_vector))
        tilt = na.arctan(r2[2]/r2[1])
        return (mag_A, mag_B, mag_C, e0_vector[0], e0_vector[1],
            e0_vector[2], tilt)


    def get_ellipsoid_parameters(self):
        r"""Calculate the parameters that describe the ellipsoid of
        the particles that constitute the halo.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        tuple : (cm, mag_A, mag_B, mag_C, e0_vector, tilt)
            The 6-tuple has in order:
              #. The center of mass as an array.
              #. mag_A as a float.
              #. mag_B as a float.
              #. mag_C as a float.
              #. e0_vector as an array.
              #. tilt as a float.
        
        Examples
        --------
        >>> params = halos[0].get_ellipsoid_parameters()
        """
        basic_parameters = self._get_ellipsoid_parameters_basic()
        toreturn = [self.center_of_mass()]
        updated = [basic_parameters[0], basic_parameters[1],
            basic_parameters[2], na.array([basic_parameters[3],
            basic_parameters[4], basic_parameters[5]]), basic_parameters[6]]
        toreturn.extend(updated)
        return tuple(toreturn)
    
    def get_ellipsoid(self):
        r"""Returns an ellipsoidal data object.
        
        This will generate a new, empty ellipsoidal data object for this
        halo.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        ellipsoid : `yt.data_objects.api.AMREllipsoidBase`
            The ellipsoidal data object.
        
        Examples
        --------
        >>> ell = halos[0].get_ellipsoid()
        """
        ep = self.get_ellipsoid_parameters()
        ell = self.data.hierarchy.ellipsoid(ep[0], ep[1], ep[2], ep[3],
            ep[4], ep[5])
        return ell
    
class HOPHalo(Halo):
    _name = "HOPHalo"
    pass


class parallelHOPHalo(Halo, ParallelAnalysisInterface):
    dont_wrap = ["maximum_density", "maximum_density_location",
        "center_of_mass", "total_mass", "bulk_velocity", "maximum_radius",
        "get_size", "get_sphere", "write_particle_list", "__getitem__",
        "virial_info", "virial_bin", "virial_mass", "virial_radius",
        "rms_velocity"]

    def virial_info(self, bins=300):
        r"""Calculates the virial information for the halo. Generally, it is
        better to call virial_radius or virial_mass instead, which calls this
        function automatically.
        """
        # Skip if we've already calculated for this number of bins.
        if self.bin_count == bins and self.overdensity is not None:
            return None
        # Do this for all because all will use it.
        self.bin_count = bins
        period = self.data.pf.domain_right_edge - \
            self.data.pf.domain_left_edge
        self.mass_bins = na.zeros(self.bin_count + 1, dtype='float64')
        cen = self.center_of_mass()
        # Cosmology
        h = self.data.pf.hubble_constant
        Om_matter = self.data.pf.omega_matter
        z = self.data.pf.current_redshift
        rho_crit = rho_crit_now * h ** 2.0 * Om_matter  # g cm^-3
        Msun2g = mass_sun_cgs
        rho_crit = rho_crit * ((1.0 + z) ** 3.0)
        # If I own some of this halo operate on the particles.
        if self.indices is not None:
            # Get some pertinent information about the halo.
            dist = na.empty(self.indices.size, dtype='float64')
            mark = 0
            # Find the distances to the particles.
            # I don't like this much, but I
            # can't see a way to eliminate a loop like this, either here or in
            # yt.math_utils.
            for pos in itertools.izip(self["particle_position_x"],
                    self["particle_position_y"], self["particle_position_z"]):
                dist[mark] = periodic_dist(cen, pos, period)
                mark += 1
            dist_min, dist_max = min(dist), max(dist)
        # If I don't have this halo, make some dummy values.
        else:
            dist_min = max(period)
            dist_max = 0.0
        # In this parallel case, we're going to find the global dist extrema
        # and built identical bins on all tasks.
        dist_min = self.comm.mpi_allreduce(dist_min, op='min')
        dist_max = self.comm.mpi_allreduce(dist_max, op='max')
        # Set up the radial bins.
        # Multiply min and max to prevent issues with digitize below.
        self.radial_bins = na.logspace(math.log10(dist_min * .99 + TINY),
            math.log10(dist_max * 1.01 + 2 * TINY), num=self.bin_count + 1)
        if self.indices is not None and self.indices.size > 1:
            # Find out which bin each particle goes into, and add the particle
            # mass to that bin.
            inds = na.digitize(dist, self.radial_bins) - 1
            for index in na.unique(inds):
                self.mass_bins[index] += \
                    na.sum(self["ParticleMassMsun"][inds == index])
            # Now forward sum the masses in the bins.
            for i in xrange(self.bin_count):
                self.mass_bins[i + 1] += self.mass_bins[i]
        # Sum up the mass_bins globally
        self.mass_bins = self.comm.mpi_allreduce(self.mass_bins, op='sum')
        # Calculate the over densities in the bins.
        self.overdensity = self.mass_bins * Msun2g / \
        (4. / 3. * math.pi * rho_crit * \
        (self.radial_bins * self.data.pf["cm"]) ** 3.0)

    def _get_ellipsoid_parameters_basic(self):
        mylog.error("Ellipsoid calculation does not work for parallelHF halos." + \
        " Please save the halos using .dump(), and reload them using" + \
        " LoadHaloes() to use this function.")
        return None

    def get_ellipsoid_parameters(self):
        r"""Calculate the parameters that describe the ellipsoid of
        the particles that constitute the halo.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        tuple : (cm, mag_A, mag_B, mag_C, e1_vector, tilt)
            The 6-tuple has in order:
              #. The center of mass as an array.
              #. mag_A as a float.
              #. mag_B as a float.
              #. mag_C as a float.
              #. e1_vector as an array.
              #. tilt as a float.
        
        Examples
        --------
        >>> params = halos[0].get_ellipsoid_parameters()
        """
        mylog.error("get_ellipsoid_parameters does not work for parallelHF halos." + \
        " Please save the halos using .dump(), and reload them using" + \
        " LoadHaloes() to use this function.")
        return None

class FOFHalo(Halo):

    def maximum_density(self):
        r"""Not implemented."""
        return -1

    def maximum_density_location(self):
        r"""Not implemented."""
        return self.center_of_mass()


class LoadedHalo(Halo):
    def __init__(self, pf, id, size=None, CoM=None,

        max_dens_point=None, group_total_mass=None, max_radius=None, bulk_vel=None,
        rms_vel=None, fnames=None, mag_A=None, mag_B=None, mag_C=None,
        e1_vec=None, tilt=None):

        self.pf = pf
        self.gridsize = (self.pf.domain_right_edge - \
            self.pf.domain_left_edge)
        self.id = id
        self.size = size
        self.CoM = CoM
        self.max_dens_point = max_dens_point
        self.group_total_mass = group_total_mass
        self.max_radius = max_radius
        self.bulk_vel = bulk_vel
        self.rms_vel = rms_vel
        self.mag_A = mag_A
        self.mag_B = mag_B
        self.mag_C = mag_C
        self.e1_vec = e1_vec
        self.tilt = tilt
        # locs=the names of the h5 files that have particle data for this halo
        self.fnames = fnames
        self.bin_count = None
        self.overdensity = None
        self.saved_fields = {}
        self.particle_mask = None
        self.ds_sort = None
        self.indices = na.array([])  # Never used for a LoadedHalo.

    def __getitem__(self, key):
        # This function will try to get particle data in one of three ways,
        # in descending preference.
        # 1. From saved_fields, e.g. we've already got it.
        # 2. From the halo h5 files off disk.
        # 3. Use the unique particle indexes of the halo to select a missing
        # field from an AMR Sphere.
        try:
            # We've already got it.
            return self.saved_fields[key]
        except KeyError:
            # Gotta go get it from the halo h5 files.
            field_data = self._get_particle_data(self.id, self.fnames,
                self.size, key)
            #if key == 'particle_position_x': field_data = None
            if field_data is not None:
                self.saved_fields[key] = field_data
                return self.saved_fields[key]
            else:
                # Dynamically create the masking array for particles, and get
                # the data using standard yt methods. The 1.05 is there to
                # account for possible silliness having to do with whether
                # the maximum density or center of mass was used to calculate
                # the maximum radius.
                ds = self.pf.h.sphere(self.CoM, 1.05 * self.max_radius)
                if self.particle_mask is None:
                    pid = self.__getitem__('particle_index')
                    sp_pid = ds['particle_index']
                    self.ds_sort = sp_pid.argsort()
                    sp_pid = sp_pid[self.ds_sort]
                    # The result of searchsorted is an array with the positions
                    # of the indexes in pid as they are in sp_pid. This is
                    # because each element of pid is in sp_pid only once.
                    self.particle_mask = na.searchsorted(sp_pid, pid)
                # We won't store this field below in saved_fields because
                # that would mean keeping two copies of it, one in the yt
                # machinery and one here.
                return ds[key][self.ds_sort][self.particle_mask]

    def _get_particle_data(self, halo, fnames, size, field):
        # Given a list of file names, a halo, its size, and the desired field,
        # this returns the particle data for that halo.
        # First get the list of fields from the first file. Not all fields
        # are saved all the time (e.g. creation_time, particle_type).
        mylog.info("Getting field %s from hdf5 halo particle files." % field)
        f = h5py.File(fnames[0], 'r')
        fields = f["Halo%08d" % halo].keys()
        # If we dont have this field, we can give up right now.
        if field not in fields:
            return None
        elif field == 'particle_index' or field == 'particle_type':
            # the only integer field
            field_data = na.empty(size, dtype='int64')
        else:
            field_data = na.empty(size, dtype='float64')
        f.close()
        # Apparently, there's a bug in h5py that was keeping the file pointer
        # f closed, even though it's re-opened below. This del seems to fix
        # that.
        del f
        offset = 0
        for fname in fnames:
            f = h5py.File(fname, 'r')
            this = f["Halo%08d" % halo][field][:]
            s = this.size
            field_data[offset:offset + s] = this
            offset += s
            f.close()
            del f
        return field_data

    def _get_ellipsoid_parameters_basic_loadedhalo(self):
        if self.mag_A is not None:
            return (self.mag_A, self.mag_B, self.mag_C, self.e1_vec[0],
                self.e1_vec[1], self.e1_vec[2], self.tilt)
        else:
            return self._get_ellipsoid_parameters_basic()

    def get_ellipsoid_parameters(self):
        r"""Calculate the parameters that describe the ellipsoid of
        the particles that constitute the halo.

        Parameters
        ----------
        None

        Returns
        -------
        tuple : (cm, mag_A, mag_B, mag_C, e1_vector, tilt)
            The 6-tuple has in order:
              #. The center of mass as an array.
              #. mag_A as a float.
              #. mag_B as a float.
              #. mag_C as a float.
              #. e1_vector as an array.
              #. tilt as a float.

        Examples
        --------
        >>> params = halos[0].get_ellipsoid_parameters()
	"""

        basic_parameters = self._get_ellipsoid_parameters_basic_loadedhalo()
        toreturn = [self.center_of_mass()]
        updated = [basic_parameters[0], basic_parameters[1],
            basic_parameters[2], na.array([basic_parameters[3],
            basic_parameters[4], basic_parameters[5]]), basic_parameters[6]]
        toreturn.extend(updated)
        return tuple(toreturn)
    
    def get_ellipsoid(self):
        r"""Returns an ellipsoidal data object.        
        This will generate a new, empty ellipsoidal data object for this
        halo.
        
        Parameters
        ----------
        None.
        
        Returns
        -------
        ellipsoid : `yt.data_objects.api.AMREllipsoidBase`
            The ellipsoidal data object.
        
        Examples
        --------
        >>> ell = halos[0].get_ellipsoid()
        """
        ep = self.get_ellipsoid_parameters()
        ell = self.pf.hierarchy.ellipsoid(ep[0], ep[1], ep[2], ep[3],
            ep[4], ep[5])
        return ell

    def get_sphere(self):
        r"""Returns a sphere source.

        This will generate a new, empty sphere source centered on this halo,
        with the maximum radius of the halo. This can be used like any other
        data container in yt.

        Parameters
        ----------
        center_of_mass : bool, optional
            True chooses the center of mass when
            calculating the maximum radius.
            False chooses from the maximum density location for HOP halos
            (it has no effect for FOF halos).
            Default = True.

        Returns
        -------
        sphere : `yt.data_objects.api.AMRSphereBase`
            The empty data source.

        Examples
        --------
        >>> sp = halos[0].get_sphere()
        """
        cen = self.center_of_mass()
        r = self.maximum_radius()
        return self.pf.h.sphere(cen, r)


class HaloList(object):

    _fields = ["particle_position_%s" % ax for ax in 'xyz']

    def __init__(self, data_source, dm_only=True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self._data_source = data_source
        self.dm_only = dm_only
        self._groups = []
        self._max_dens = {}
        self.__obtain_particles()
        self._run_finder()
        mylog.info("Parsing outputs")
        self._parse_output()
        mylog.debug("Finished. (%s)", len(self))

    def __obtain_particles(self):
        if self.dm_only:
            ii = self._get_dm_indices()
        else:
            ii = slice(None)
        self.particle_fields = {}
        for field in self._fields:
            tot_part = self._data_source[field].size
            if field == "particle_index":
                self.particle_fields[field] = \
                    self._data_source[field][ii].astype('int64')
            else:
                self.particle_fields[field] = \
                    self._data_source[field][ii].astype('float64')
            del self._data_source[field]
        self._base_indices = na.arange(tot_part)[ii]
        gc.collect()

    def _get_dm_indices(self):
        if 'creation_time' in self._data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on creation time")
            return (self._data_source["creation_time"] < 0)
        elif 'particle_type' in self._data_source.hierarchy.field_list:
            mylog.debug("Differentiating based on particle type")
            return (self._data_source["particle_type"] == 1)
        else:
            mylog.warning("No particle_type, no creation_time, so not distinguishing.")
            return slice(None)

    def _parse_output(self):
        unique_ids = na.unique(self.tags)
        counts = na.bincount(self.tags + 1)
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        dens = self.densities[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i + 1]
            if i == -1:
                cp += counts[i + 1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(self._halo_class(self, i, group_indices))
            md_i = na.argmax(dens[cp:cp_c])
            px, py, pz = \
                [self.particle_fields['particle_position_%s' % ax][group_indices]
                                            for ax in 'xyz']
            self._max_dens[i] = (dens[cp:cp_c][md_i], px[md_i],
                py[md_i], pz[md_i])
            cp += counts[i + 1]

    def __len__(self):
        return len(self._groups)

    def __iter__(self):
        for i in self._groups:
            yield i

    def __getitem__(self, key):
        return self._groups[key]

    def nearest_neighbors_3D(self, haloID, num_neighbors=7, search_radius=.2):
        r"""For a halo its nearest neighbors in 3D using the kd tree.

        This will calculate the nearest neighbors of a halo, using the kD tree.
        Returns a list of the neighbors distances and ID with format
        [distance,haloID].

        Parameters
        ----------
        haloID : integer
            The halo to find neighbors for.
        num_neighbors : integer
            How many neighbors to search for. Default = 7.
        search_radius : float
            How far away to look for neighbors in code units. Default = 0.2.

        Examples
        --------
        >>> neighbors = halos.nearest_neighbors_3D(0)
        """
        period = self.pf.domain_right_edge - self.pf.domain_left_edge
        # Initialize the dataset of points from all the haloes
        dataset = []
        for group in self:
            p = Point()
            p.data = group.center_of_mass().tolist()
            p.haloID = group.id
            dataset.append(p)
        mylog.info('Building kd tree...')
        kd = buildKdHyperRectTree(dataset[:], 2 * num_neighbors)
        # make the neighbors object
        neighbors = Neighbors()
        neighbors.k = num_neighbors
        neighbors.points = []
        neighbors.minDistanceSquared = search_radius * search_radius
        mylog.info('Finding nearest neighbors...')
        getKNN(self[haloID].center_of_mass().tolist(), kd, neighbors, 0.,
            period.tolist())
        # convert the data in order to return something less perverse than a
        # Neighbors object, also root the distances
        n_points = []
        for n in neighbors.points:
            n_points.append([math.sqrt(n[0]), n[1].haloID])
        return n_points

    def nearest_neighbors_2D(self, haloID, num_neighbors=7, search_radius=.2,
        proj_dim=0):
        r"""For a halo its nearest neighbors in 2D using the kd tree.

        This will strip a dimension from consideration in the kD-tree, and then
        calculate all the nearest projected neighbors of a halo.  Returns a
        list of the neighbors distances and ID with format [distance,haloID].

        Parameters
        ----------
        haloID : int
            The halo to find neighbors for.
        num_neighbors : int
            How many neighbors to search for. Default = 7.
        search_radius : float
            How far away to look for neighbors in code units. Default = 0.2.
        proj_dim : int
            Which dimension (0, 1, or 2) to project the halos into 2D.
            Default = 0.

        Examples
        --------
        >>> neighbors = halos.nearest_neighbors_2D(0)
        """
        # Set up a vector to multiply other
        # vectors by to project along proj_dim
        vec = na.array([1., 1., 1.])
        vec[proj_dim] = 0.
        period = self.pf.domain_right_edge - self.pf.domain_left_edge
        period = period * vec
        # Initialize the dataset of points from all the haloes
        dataset = []
        for group in self:
            p = Point()
            cm = group.center_of_mass() * vec
            p.data = cm.tolist()
            p.haloID = group.id
            dataset.append(p)
        mylog.info('Building kd tree...')
        kd = buildKdHyperRectTree(dataset[:], 2 * num_neighbors)
        # make the neighbors object
        neighbors = Neighbors()
        neighbors.k = num_neighbors
        neighbors.points = []
        neighbors.minDistanceSquared = search_radius * search_radius
        mylog.info('Finding nearest neighbors...')
        cm = self[haloID].center_of_mass() * vec
        getKNN(cm.tolist(), kd, neighbors, 0., period.tolist())
        # convert the data in order to return something less perverse than a
        # Neighbors object, also root the distances
        n_points = []
        for n in neighbors.points:
            n_points.append([math.sqrt(n[0]), n[1].haloID])
        return n_points

    def write_out(self, filename, ellipsoid_data=False):
        r"""Write out standard halo information to a text file.

        Parameters
        ----------
        filename : String
            The name of the file to write to.

        ellipsoid_data : bool.
            Whether to print the ellipsoidal information to the file.
            Default = False.
        
        Examples
        --------
        >>> halos.write_out("HopAnalysis.out")
        """
        if hasattr(filename, 'write'):
            f = filename
        else:
            f = open(filename, "w")
        f.write("# HALOS FOUND WITH %s\n" % (self._name))

        if not ellipsoid_data:
            f.write("\t".join(["# Group","Mass","# part","max dens"
                               "x","y","z", "center-of-mass",
                               "x","y","z",
                               "vx","vy","vz","max_r","rms_v","\n"]))
        else:
            f.write("\t".join(["# Group","Mass","# part","max dens"
                               "x","y","z", "center-of-mass",
                               "x","y","z",
                               "vx","vy","vz","max_r","rms_v",
                               "mag_A", "mag_B", "mag_C", "e1_vec0",
                               "e1_vec1", "e1_vec2", "tilt", "\n"]))

        for group in self:
            f.write("%10i\t" % group.id)
            f.write("%0.9e\t" % group.total_mass())
            f.write("%10i\t" % group.get_size())
            f.write("%0.9e\t" % group.maximum_density())
            f.write("\t".join(["%0.9e" % v for v in \
                group.maximum_density_location()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.center_of_mass()]))
            f.write("\t")
            f.write("\t".join(["%0.9e" % v for v in group.bulk_velocity()]))
            f.write("\t")
            f.write("%0.9e\t" % group.maximum_radius())
            f.write("%0.9e\t" % group.rms_velocity())
            if ellipsoid_data:
                f.write("\t".join(["%0.9e" % v for v in group._get_ellipsoid_parameters_basic()]))
            f.write("\n")
            f.flush()
        f.close()

    def write_particle_lists_txt(self, prefix, fp=None):
        r"""Write out the names of the HDF5 files containing halo particle data
        to a text file. Needed in particular for parallel analysis output.

        Parameters
        ----------
        prefix : String
            The prefix for the name of the file.

        Examples
        --------
        >>> halos.write_particle_lists_txt("halo-parts")
        """
        if hasattr(fp, 'write'):
            f = fp
        else:
            f = open("%s.txt" % prefix, "w")
        for group in self:
            if group.tasks is not None:
                fn = ""
                for task in group.tasks:
                    fn += "%s.h5 " % self.comm.get_filename(prefix, rank=task)
            elif self._distributed:
                fn = "%s.h5" % self.comm.get_filename(prefix,
                    rank=group._owner)
            else:
                fn = "%s.h5" % self.comm.get_filename(prefix)
            gn = "Halo%08i" % (group.id)
            f.write("%s %s\n" % (gn, fn))
            f.flush()
        f.close()


class HOPHaloList(HaloList):

    _name = "HOP"
    _halo_class = HOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["ParticleMassMsun"]

    def __init__(self, data_source, threshold=160.0, dm_only=True):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        self.threshold = threshold
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        self.densities, self.tags = \
            RunHOP(self.particle_fields["particle_position_x"] / self.period[0],
                self.particle_fields["particle_position_y"] / self.period[1],
                self.particle_fields["particle_position_z"] / self.period[2],
                self.particle_fields["ParticleMassMsun"],
                self.threshold)
        self.particle_fields["densities"] = self.densities
        self.particle_fields["tags"] = self.tags

    def write_out(self, filename="HopAnalysis.out", ellipsoid_data=False):
        r"""Write out standard halo information to a text file.

        Parameters
        ----------
        filename : String
            The name of the file to write to. Default = "HopAnalysis.out".

        ellipsoid_data : bool.
            Whether to print the ellipsoidal information to the file.
            Default = False.

        Examples
        --------
        >>> halos.write_out("HopAnalysis.out")
        """
        HaloList.write_out(self, filename, ellipsoid_data)


class FOFHaloList(HaloList):
    _name = "FOF"
    _halo_class = FOFHalo

    def __init__(self, data_source, link=0.2, dm_only=True):
        self.link = link
        mylog.info("Initializing FOF")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        self.tags = \
        RunFOF(self.particle_fields["particle_position_x"] / self.period[0],
               self.particle_fields["particle_position_y"] / self.period[1],
               self.particle_fields["particle_position_z"] / self.period[2],
               self.link)
        self.densities = na.ones(self.tags.size, dtype='float64') * -1
        self.particle_fields["densities"] = self.densities
        self.particle_fields["tags"] = self.tags

    def write_out(self, filename="FOFAnalysis.out", ellipsoid_data=False):
        r"""Write out standard halo information to a text file.

        Parameters
        ----------
        filename : String
            The name of the file to write to. Default = "FOFAnalysis.out".

        ellipsoid_data : bool.
            Whether to print the ellipsoidal information to the file.
            Default = False.

        Examples
        --------
        >>> halos.write_out("FOFAnalysis.out")
        """
        HaloList.write_out(self, filename, ellipsoid_data)


class LoadedHaloList(HaloList):
    _name = "Loaded"

    def __init__(self, pf, basename):
        ParallelAnalysisInterface.__init__(self)
        self.pf = pf
        self._groups = []
        self.basename = basename
        self._retrieve_halos()

    def _retrieve_halos(self):
        # First get the halo particulars.
        lines = file("%s.out" % self.basename)
        # The location of particle data for each halo.
        locations = self._collect_halo_data_locations()
        halo = 0
        for line in lines:
            orig = line
            # Skip the comment lines at top.
            if line[0] == "#": continue
            line = line.split()
            # get the particle data
            size = int(line[2])
            fnames = locations[halo]
            # Everything else
            CoM = na.array([float(line[7]), float(line[8]), float(line[9])])
            max_dens_point = na.array([float(line[3]), float(line[4]),
                float(line[5]), float(line[6])])
            group_total_mass = float(line[1])
            max_radius = float(line[13])
            bulk_vel = na.array([float(line[10]), float(line[11]),
                float(line[12])])
            rms_vel = float(line[14])
            if len(line) == 15:
                # No ellipsoid information
                self._groups.append(LoadedHalo(self.pf, halo, size = size,
                    CoM = CoM,
                    max_dens_point = max_dens_point,
                    group_total_mass = group_total_mass, max_radius = max_radius,
                    bulk_vel = bulk_vel, rms_vel = rms_vel, fnames = fnames))
            elif len(line) == 22:
                # Ellipsoid information
                mag_A = float(line[15])
                mag_B = float(line[16])
                mag_C = float(line[17])
                e1_vec0 = float(line[18])
                e1_vec1 = float(line[19])
                e1_vec2 = float(line[20])
                e1_vec = na.array([e1_vec0, e1_vec1, e1_vec2])
                tilt = float(line[21])
                self._groups.append(LoadedHalo(self.pf, halo, size = size,
                    CoM = CoM,
                    max_dens_point = max_dens_point,
                    group_total_mass = group_total_mass, max_radius = max_radius,
                    bulk_vel = bulk_vel, rms_vel = rms_vel, fnames = fnames,
                    mag_A = mag_A, mag_B = mag_B, mag_C = mag_C, e1_vec = e1_vec,
                    tilt = tilt))
            else:
                mylog.error("I am unable to parse this line. Too many or too few items. %s" % orig)
            halo += 1

    def _collect_halo_data_locations(self):
        # The halos are listed in order in the file.
        lines = file("%s.txt" % self.basename)
        locations = []
        realpath = path.realpath("%s.txt" % self.basename)
        for line in lines:
            line = line.split()
            # Prepend the hdf5 file names with the full path.
            temp = []
            for item in line[1:]:
                # This assumes that the .txt is in the same place as
                # the h5 files, which is a good one I think.
                item = item.split("/")
                temp.append(path.join(path.dirname(realpath), item[-1]))
            locations.append(temp)
        lines.close()
        return locations


class parallelHOPHaloList(HaloList, ParallelAnalysisInterface):
    _name = "parallelHOP"
    _halo_class = parallelHOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["ParticleMassMsun", "particle_index"]

    def __init__(self, data_source, padding, num_neighbors, bounds, total_mass,
        period, threshold=160.0, dm_only=True, rearrange=True, premerge=True,
        tree='F'):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is set, only run it on the dark matter particles, otherwise
        on all particles.  Returns an iterable collection of *HopGroup* items.
        """
        ParallelAnalysisInterface.__init__(self)
        self.threshold = threshold
        self.num_neighbors = num_neighbors
        self.bounds = bounds
        self.total_mass = total_mass
        self.rearrange = rearrange
        self.period = period
        self.old_period = period.copy()
        self.period = na.array([1.] * 3)
        self._data_source = data_source
        self.premerge = premerge
        self.tree = tree
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only)

    def _run_finder(self):
        yt_counters("Reading Data")
        # Test to make sure the particle IDs aren't suspicious.
        exit = False
        if (self.particle_fields["particle_index"] < 0).any():
            mylog.error("Negative values in particle_index field. Parallel HOP will fail.")
            exit = True
        if na.unique(self.particle_fields["particle_index"]).size != \
                self.particle_fields["particle_index"].size:
            mylog.error("Non-unique values in particle_index field. Parallel HOP will fail.")
            exit = True

        self.comm.mpi_exit_test(exit)
        # Try to do this in a memory conservative way.
        na.divide(self.particle_fields['ParticleMassMsun'], self.total_mass,
            self.particle_fields['ParticleMassMsun'])
        na.divide(self.particle_fields["particle_position_x"],
            self.old_period[0], self.particle_fields["particle_position_x"])
        na.divide(self.particle_fields["particle_position_y"],
            self.old_period[1], self.particle_fields["particle_position_y"])
        na.divide(self.particle_fields["particle_position_z"],
            self.old_period[2], self.particle_fields["particle_position_z"])
        obj = ParallelHOPHaloFinder(self.period, self.padding,
            self.num_neighbors, self.bounds,
            self.particle_fields,
            self.threshold, rearrange=self.rearrange, premerge=self.premerge,
            tree=self.tree)
        self.densities, self.tags = obj.density, obj.chainID
        # I'm going to go ahead and delete self.densities because it's not
        # actually being used. I'm not going to remove it altogether because
        # it may be useful to someone someday.
        del self.densities
        self.group_count = obj.group_count
        self.group_sizes = obj.group_sizes
        if self.group_count == 0:
            mylog.info("There are no halos found.")
            return
        self.CoM = obj.CoM
        self.Tot_M = obj.Tot_M * self.total_mass
        self.max_dens_point = obj.max_dens_point
        self.max_radius = obj.max_radius
        for dd in range(3):
            self.CoM[:, dd] *= self.old_period[dd]
            self.max_dens_point[:, dd + 1] *= self.old_period[dd]
        # This is wrong, below, with uneven boundaries. We'll cross that bridge
        # when we get there.
        self.max_radius *= self.old_period[0]
        self.period = self.old_period.copy()
        # Precompute the bulk velocity in parallel.
        yt_counters("Precomp bulk vel.")
        self.bulk_vel = na.zeros((self.group_count, 3), dtype='float64')
        yt_counters("bulk vel. reading data")
        pm = obj.mass
        # Fix this back to un-normalized units.
        na.multiply(pm, self.total_mass, pm)
        xv = self._data_source["particle_velocity_x"][self._base_indices]
        yv = self._data_source["particle_velocity_y"][self._base_indices]
        zv = self._data_source["particle_velocity_z"][self._base_indices]
        yt_counters("bulk vel. reading data")
        yt_counters("bulk vel. computing")
        select = (self.tags >= 0)
        calc = len(na.where(select == True)[0])
        if calc:
            vel = na.empty((calc, 3), dtype='float64')
            ms = pm[select]
            vel[:, 0] = xv[select] * ms
            vel[:, 1] = yv[select] * ms
            vel[:, 2] = zv[select] * ms
            subchain = self.tags[select]
            sort = subchain.argsort()
            vel = vel[sort]
            sort_subchain = subchain[sort]
            uniq_subchain = na.unique(sort_subchain)
            diff_subchain = na.ediff1d(sort_subchain)
            marks = (diff_subchain > 0)
            marks = na.arange(calc)[marks] + 1
            marks = na.concatenate(([0], marks, [calc]))
            for i, u in enumerate(uniq_subchain):
                self.bulk_vel[u] = na.sum(vel[marks[i]:marks[i + 1]], axis=0)
            del vel, subchain, sort_subchain
            del diff_subchain
        # Bring it together, and divide by the previously computed total mass
        # of each halo.
        self.bulk_vel = self.comm.mpi_allreduce(self.bulk_vel, op='sum')
        for groupID in xrange(self.group_count):
            self.bulk_vel[groupID] = \
                self.bulk_vel[groupID] / self.Tot_M[groupID]
        yt_counters("bulk vel. computing")
        # Now calculate the RMS velocity of the groups in parallel, very
        # similarly to the bulk velocity and re-using some of the arrays.
        yt_counters("rms vel computing")
        rms_vel_temp = na.zeros((self.group_count, 2), dtype='float64')
        if calc:
            vel = na.empty((calc, 3), dtype='float64')
            vel[:, 0] = xv[select] * ms
            vel[:, 1] = yv[select] * ms
            vel[:, 2] = zv[select] * ms
            vel = vel[sort]
            for i, u in enumerate(uniq_subchain):
                # This finds the sum locally.
                rms_vel_temp[u][0] = na.sum(((vel[marks[i]:marks[i + 1]] - \
                    self.bulk_vel[u]) / self.Tot_M[u]) ** 2.)
                # I could use self.group_sizes...
                rms_vel_temp[u][1] = marks[i + 1] - marks[i]
            del vel, marks, uniq_subchain
        # Bring it together.
        rms_vel_temp = self.comm.mpi_allreduce(rms_vel_temp, op='sum')
        self.rms_vel = na.empty(self.group_count, dtype='float64')
        for groupID in xrange(self.group_count):
            # Here we do the Mean and the Root.
            self.rms_vel[groupID] = \
                na.sqrt(rms_vel_temp[groupID][0] / rms_vel_temp[groupID][1]) * \
                self.group_sizes[groupID]
        del rms_vel_temp
        yt_counters("rms vel computing")
        self.taskID = obj.mine
        self.halo_taskmap = obj.halo_taskmap  # A defaultdict.
        del obj
        gc.collect()
        yt_counters("Precomp bulk vel.")

    def _parse_output(self):
        yt_counters("Final Grouping")
        """
        Each task will make an entry for all groups, but it may be empty.
        """
        unique_ids = na.unique(self.tags)
        counts = na.bincount((self.tags + 1).tolist())
        sort_indices = na.argsort(self.tags)
        grab_indices = na.indices(self.tags.shape).ravel()[sort_indices]
        del sort_indices
        cp = 0
        index = 0
        # We want arrays for parallel HOP
        self._groups = na.empty(self.group_count, dtype='object')
        self._max_dens = na.empty((self.group_count, 4), dtype='float64')
        if self.group_count == 0:
            mylog.info("There are no halos found.")
            return
        for i in unique_ids:
            if i == -1:
                cp += counts[i + 1]
                continue
            # If there is a gap in the unique_ids, make empty groups to
            # fill it in.
            while index < i:
                self._groups[index] = self._halo_class(self, index, \
                    size=self.group_sizes[index], CoM=self.CoM[index], \
                    max_dens_point=self.max_dens_point[index], \
                    group_total_mass=self.Tot_M[index],
                    max_radius=self.max_radius[index],
                    bulk_vel=self.bulk_vel[index],
                    tasks=self.halo_taskmap[index],
                    rms_vel=self.rms_vel[index])
                # I don't own this halo
                self.comm.do_not_claim_object(self._groups[index])
                self._max_dens[index] = [self.max_dens_point[index][0],
                    self.max_dens_point[index][1], \
                    self.max_dens_point[index][2],
                    self.max_dens_point[index][3]]
                index += 1
            cp_c = cp + counts[i + 1]
            group_indices = grab_indices[cp:cp_c]
            self._groups[index] = self._halo_class(self, i, group_indices, \
                size=self.group_sizes[i], CoM=self.CoM[i], \
                max_dens_point=self.max_dens_point[i], \
                group_total_mass=self.Tot_M[i], max_radius=self.max_radius[i],
                bulk_vel=self.bulk_vel[i], tasks=self.halo_taskmap[index],
                rms_vel=self.rms_vel[i])
            # This halo may be owned by many, including this task
            self.comm.claim_object(self._groups[index])
            self._max_dens[index] = [self.max_dens_point[i][0],
                self.max_dens_point[i][1], \
                self.max_dens_point[i][2], self.max_dens_point[i][3]]
            cp += counts[i + 1]
            index += 1
        # If there are missing groups at the end, add them.
        while index < self.group_count:
            self._groups[index] = self._halo_class(self, index, \
                size=self.group_sizes[index], CoM=self.CoM[index], \
                max_dens_point=self.max_dens_point[i], \
                group_total_mass=self.Tot_M[index],
                max_radius=self.max_radius[index],
                bulk_vel=self.bulk_vel[index], tasks=self.halo_taskmap[index],
                rms_vel=self.rms_vel[index])
            self.comm.do_not_claim_object(self._groups[index])
            self._max_dens[index] = [self.max_dens_point[index][0],
                self.max_dens_point[index][1], \
                self.max_dens_point[index][2], self.max_dens_point[index][3]]
            index += 1
        # Clean up
        del self.max_dens_point, self.max_radius, self.bulk_vel
        del self.halo_taskmap, self.tags, self.rms_vel
        del grab_indices, unique_ids, counts
        try:
            del group_indices
        except UnboundLocalError:
            pass

    def __len__(self):
        return self.group_count

    def write_out(self, filename="parallelHopAnalysis.out", ellipsoid_data=False):
        r"""Write out standard halo information to a text file.

        Parameters
        ----------
        filename : String
            The name of the file to write to.
            Default = "parallelHopAnalysis.out".

        Examples
        --------
        >>> halos.write_out("parallelHopAnalysis.out")
        """
        HaloList.write_out(self, filename, ellipsoid_data)


class GenericHaloFinder(HaloList, ParallelAnalysisInterface):
    def __init__(self, pf, ds, dm_only=True, padding=0.0):
        ParallelAnalysisInterface.__init__(self)
        self.pf = pf
        self.hierarchy = pf.h
        self.center = (na.array(ds.right_edge) + na.array(ds.left_edge)) / 2.0

    def _parse_halolist(self, threshold_adjustment):
        groups = []
        max_dens = {}
        hi = 0
        LE, RE = self.bounds
        for halo in self._groups:
            this_max_dens = halo.maximum_density_location()
            # if the most dense particle is in the box, keep it
            if na.all((this_max_dens >= LE) & (this_max_dens <= RE)):
                # Now we add the halo information to OURSELVES, taken from the
                # self.hop_list
                # We need to mock up the HOPHaloList thingie, so we need to
                #     set self._max_dens
                max_dens_temp = list(self._max_dens[halo.id])[0] / \
                    threshold_adjustment
                max_dens[hi] = [max_dens_temp] + \
                    list(self._max_dens[halo.id])[1:4]
                groups.append(self._halo_class(self, hi))
                groups[-1].indices = halo.indices
                self.comm.claim_object(groups[-1])
                hi += 1
        del self._groups, self._max_dens  # explicit >> implicit
        self._groups = groups
        self._max_dens = max_dens

    def _join_halolists(self):
        # First we get the total number of halos the entire collection
        # has identified
        # Note I have added a new method here to help us get information
        # about processors and ownership and so forth.
        # _mpi_info_dict returns a dict of {proc: whatever} where whatever is
        # what is fed in on each proc.
        mine, halo_info = self.comm.mpi_info_dict(len(self))
        nhalos = sum(halo_info.values())
        # Figure out our offset
        my_first_id = sum([v for k, v in halo_info.items() if k < mine])
        # Fix our max_dens
        max_dens = {}
        for i, m in self._max_dens.items():
            max_dens[i + my_first_id] = m
        self._max_dens = max_dens
        for halo in self._groups:
            halo._max_dens = self._max_dens
        # sort the list by the size of the groups
        # Now we add ghost halos and reassign all the IDs
        # Note: we already know which halos we own!
        after = my_first_id + len(self._groups)
        # One single fake halo, not owned, does the trick
        self._groups = [self._halo_class(self, i) for i in range(my_first_id)] + \
                       self._groups + \
                       [self._halo_class(self, i) for i in range(after, nhalos)]
        id = 0
        for proc in sorted(halo_info.keys()):
            for halo in self._groups[id:id + halo_info[proc]]:
                halo.id = id
                halo._distributed = self._distributed
                halo._owner = proc
                id += 1

        def haloCmp(h1, h2):
            c = cmp(h1.total_mass(), h2.total_mass())
            if c != 0:
                return -1 * c
            if c == 0:
                return cmp(h1.center_of_mass()[0], h2.center_of_mass()[0])
        self._groups.sort(haloCmp)
        sorted_max_dens = {}
        for i, halo in enumerate(self._groups):
            if halo.id in self._max_dens:
                sorted_max_dens[i] = self._max_dens[halo.id]
            halo.id = i
        self._max_dens = sorted_max_dens
        for i, halo in enumerate(self._groups):
            halo._max_dens = self._max_dens

    def _reposition_particles(self, bounds):
        # This only does periodicity.  We do NOT want to deal with anything
        # else.  The only reason we even do periodicity is the
        LE, RE = bounds
        dw = self.pf.domain_right_edge - self.pf.domain_left_edge
        for i, ax in enumerate('xyz'):
            arr = self._data_source["particle_position_%s" % ax]
            arr[arr < LE[i] - self.padding] += dw[i]
            arr[arr > RE[i] + self.padding] -= dw[i]

    def write_out(self, filename, ellipsoid_data=False):
        r"""Write out standard halo information to a text file.

        Parameters
        ----------
        filename : String
            The name of the file to write to.

        ellipsoid_data : bool.
            Whether to print the ellipsoidal information to the file.
            Default = False.

        Examples
        --------
        >>> halos.write_out("HopAnalysis.out")
        """
        f = self.comm.write_on_root(filename)
        HaloList.write_out(self, f, ellipsoid_data)

    def write_particle_lists_txt(self, prefix):
        r"""Write out the names of the HDF5 files containing halo particle data
        to a text file.

        This function wirtes out the names of all the HDF5 files that would
        contain halo particle data.  Only the root processor writes out.

        Parameters
        ----------
        prefix : String
            The prefix for the name of the file.

        Examples
        --------
        >>> halos.write_particle_lists_txt("halo-parts")
        """
        f = self.comm.write_on_root("%s.txt" % prefix)
        HaloList.write_particle_lists_txt(self, prefix, fp=f)

    @parallel_blocking_call
    def write_particle_lists(self, prefix):
        r"""Write out the particle data for halos to HDF5 files.

        This function will accept a filename prefix, and for every halo it will
        write out an HDF5 file containing the positions, velocities, indices
        and masses of the constituent particles.  However, if the halo finder
        is run in parallel, halos will only be written out on the processors to
        which they belong.  See `Halo.write_particle_lists_txt` for how to
        track these halos globally across files.

        Parameters
        ----------
        prefix : String
            The prefix for the name(s) of the HDF5 files.

        Examples
        --------
        >>> halos.write_particle_lists("halo-parts")
        """
        fn = "%s.h5" % self.comm.get_filename(prefix)
        f = h5py.File(fn, "w")
        for halo in self._groups:
            if not self.comm.is_mine(halo): continue
            halo.write_particle_list(f)
        f.close()

    def dump(self, basename="HopAnalysis", ellipsoid_data=False):
        r"""Save the full halo data to disk.

        This function will save the halo data in such a manner that it can be
        easily re-loaded later using `GenericHaloFinder.load`.
        This is similar in concept to
        pickling the data, but outputs the data in the already-established
        data formats. The simple halo data is written to a text file
        (e.g. "HopAnalysis.out") using write_out(), and the particle data
        to hdf5 files (e.g. "HopAnalysis.h5")
        using write_particle_lists().

        Parameters
        ----------
        basename : String
            The base name for the files the data will be written to. Default =
            "HopAnalysis".

        ellipsoid_data : bool.
            Whether to save the ellipsoidal information to the files.
            Default = False.
        
        Examples
        --------
        >>> halos.dump("MyHalos")
        """
        self.write_out("%s.out" % basename, ellipsoid_data)
        self.write_particle_lists(basename)
        self.write_particle_lists_txt(basename)


class parallelHF(GenericHaloFinder, parallelHOPHaloList):
    def __init__(self, pf, subvolume=None, threshold=160, dm_only=True, \
        resize=True, rearrange=True,\
        fancy_padding=True, safety=1.5, premerge=True, sample=0.03, \
        total_mass=None, num_particles=None, tree='F'):
        r"""Parallel HOP halo finder.

        Halos are built by:
        1. Calculating a density for each particle based on a smoothing kernel.
        2. Recursively linking particles to other particles from lower density
        particles to higher.
        3. Geometrically proximate chains are identified and
        4. merged into final halos following merging rules.

        Lower thresholds generally produce more halos, and the largest halos
        become larger. Also, halos become more filamentary and over-connected.

        This is very similar to HOP, but it does not produce precisely the
        same halos due to unavoidable numerical differences.

        Skory et al. "Parallel HOP: A Scalable Halo Finder for Massive
        Cosmological Data Sets." arXiv (2010) 1001.3411

        Parameters
        ----------
        pf : `StaticOutput`
            The parameter file on which halo finding will be conducted.
        threshold : float
            The density threshold used when building halos. Default = 160.0.
        dm_only : bool
            If True, only dark matter particles are used when building halos.
            Default = False.
        resize : bool
            Turns load-balancing on or off. Default = True.
        kdtree : string
            Chooses which kD Tree to use. The Fortran one (kdtree = 'F') is
            faster, but uses more memory. The Cython one (kdtree = 'C') is
            slower but is more memory efficient.
            Default = 'F'
        rearrange : bool
            Turns on faster nearest neighbor searches at the cost of increased
            memory usage.
            This option only applies when using the Fortran tree.
            Default = True.
        fancy_padding : bool
            True calculates padding independently for each face of each
            subvolume. Default = True.
        safety : float
            Due to variances in inter-particle spacing in the volume, the
            padding may need to be increased above the raw calculation.
            This number is multiplied to the calculated padding, and values
            >1 increase the padding. Default = 1.5.
        premerge : bool
            True merges chains in two steps (rather than one with False), which
            can speed up halo finding by 25% or more. However, True can result
            in small (<<1%) variations in the final halo masses when compared
            to False. Default = True.
        sample : float
            The fraction of the full dataset on which load-balancing is
            performed. Default = 0.03.
        total_mass : float
            If HOP is run on the same dataset mulitple times, the total mass
            of particles in Msun units in the full volume can be supplied here
            to save time.
            This must correspond to the particles being operated on, meaning
            if stars are included in the halo finding, they must be included
            in this mass as well, and visa-versa.
            If halo finding on a subvolume, this still corresponds with the
            mass in the entire volume.
            Default = None, which means the total mass is automatically
            calculated.
        num_particles : integer
            The total number of particles in the volume, in the same fashion
            as `total_mass` is calculated. Specifying this turns off
            fancy_padding.
            Default = None, which means the number of particles is
            automatically calculated.

        Examples
        -------
        >>> pf = load("RedshiftOutput0000")
        >>> halos = parallelHF(pf)
        """
        if subvolume is not None:
            ds_LE = na.array(subvolume.left_edge)
            ds_RE = na.array(subvolume.right_edge)
        self._data_source = pf.h.all_data()
        GenericHaloFinder.__init__(self, pf, self._data_source, dm_only,
            padding=0.0)
        self.padding = 0.0
        self.num_neighbors = 65
        self.safety = safety
        self.sample = sample
        self.tree = tree
        if self.tree != 'F' and self.tree != 'C':
            mylog.error("No kD Tree specified!")
        period = pf.domain_right_edge - pf.domain_left_edge
        topbounds = na.array([[0., 0., 0.], period])
        # Cut up the volume evenly initially, with no padding.
        padded, LE, RE, self._data_source = \
            self.partition_hierarchy_3d(ds=self._data_source,
            padding=self.padding)
        # also get the total mass of particles
        yt_counters("Reading Data")
        # Adaptive subregions by bisection. We do not load balance if we are
        # analyzing a subvolume.
        ds_names = ["particle_position_x", "particle_position_y",
            "particle_position_z"]
        if ytcfg.getboolean("yt", "inline") == False and \
            resize and self.comm.size != 1 and subvolume is None:
            random.seed(self.comm.rank)
            cut_list = self.partition_hierarchy_3d_bisection_list()
            root_points = self._subsample_points()
            self.bucket_bounds = []
            if self.comm.rank == 0:
                self._recursive_divide(root_points, topbounds, 0, cut_list)
            self.bucket_bounds = \
                self.comm.mpi_bcast(self.bucket_bounds)
            my_bounds = self.bucket_bounds[self.comm.rank]
            LE, RE = my_bounds[0], my_bounds[1]
            self._data_source = self.hierarchy.region_strict([0.] * 3, LE, RE)
        # If this isn't parallel, define the region as an AMRRegionStrict so
        # particle IO works.
        if self.comm.size == 1:
            self._data_source = self.hierarchy.periodic_region_strict([0.5] * 3,
                LE, RE)
        # get the average spacing between particles for this region
        # The except is for the serial case where the full box is what we want.
        if num_particles is None:
            data = self._data_source["particle_position_x"]
        try:
            l = self._data_source.right_edge - self._data_source.left_edge
        except AttributeError:
            l = pf.domain_right_edge - pf.domain_left_edge
        vol = l[0] * l[1] * l[2]
        full_vol = vol
        # We will use symmetric padding when a subvolume is being used.
        if not fancy_padding or subvolume is not None or \
                num_particles is not None:
            if num_particles is None:
                num_particles = data.size
            avg_spacing = (float(vol) / num_particles) ** (1. / 3.)
            # padding is a function of inter-particle spacing, this is an
            # approximation, but it's OK with the safety factor
            padding = (self.num_neighbors) ** (1. / 3.) * self.safety * \
                avg_spacing
            self.padding = (na.ones(3, dtype='float64') * padding,
                na.ones(3, dtype='float64') * padding)
            mylog.info('padding %s avg_spacing %f vol %f local_parts %d' % \
                (str(self.padding), avg_spacing, vol, num_particles))
        # Another approach to padding, perhaps more accurate.
        elif fancy_padding and self._distributed:
            LE_padding = na.empty(3, dtype='float64')
            RE_padding = na.empty(3, dtype='float64')
            avg_spacing = (float(vol) / data.size) ** (1. / 3.)
            base_padding = (self.num_neighbors) ** (1. / 3.) * self.safety * \
                avg_spacing
            for dim in xrange(3):
                if ytcfg.getboolean("yt", "inline") == False:
                    data = self._data_source[ds_names[dim]]
                else:
                    data = self._data_source[ds_names[dim]]
                num_bins = 1000
                width = self._data_source.right_edge[dim] - \
                    self._data_source.left_edge[dim]
                area = (self._data_source.right_edge[(dim + 1) % 3] - \
                    self._data_source.left_edge[(dim + 1) % 3]) * \
                    (self._data_source.right_edge[(dim + 2) % 3] - \
                    self._data_source.left_edge[(dim + 2) % 3])
                bin_width = base_padding
                num_bins = int(math.ceil(width / bin_width))
                bins = na.arange(num_bins + 1, dtype='float64') * bin_width + \
                    self._data_source.left_edge[dim]
                counts, bins = na.histogram(data, bins)
                # left side.
                start = 0
                count = counts[0]
                while count < self.num_neighbors:
                    start += 1
                    count += counts[start]
                # Get the avg spacing in just this boundary.
                vol = area * (bins[start + 1] - bins[0])
                avg_spacing = (float(vol) / count) ** (1. / 3.)
                LE_padding[dim] = (self.num_neighbors) ** (1. / 3.) * \
                    self.safety * avg_spacing
                # right side.
                start = -1
                count = counts[-1]
                while count < self.num_neighbors:
                    start -= 1
                    count += counts[start]
                vol = area * (bins[-1] - bins[start - 1])
                avg_spacing = (float(vol) / count) ** (1. / 3.)
                RE_padding[dim] = (self.num_neighbors) ** (1. / 3.) * \
                    self.safety * avg_spacing
            self.padding = (LE_padding, RE_padding)
            del bins, counts
            mylog.info('fancy_padding %s avg_spacing %f full_vol %f local_parts %d %s' % \
                (str(self.padding), avg_spacing, full_vol,
                data.size, str(self._data_source)))
        # Now we get the full box mass after we have the final composition of
        # subvolumes.
        if total_mass is None:
            total_mass = self.comm.mpi_allreduce((self._data_source["ParticleMassMsun"].astype('float64')).sum(),
                                                 op='sum')
        if not self._distributed:
            self.padding = (na.zeros(3, dtype='float64'),
                na.zeros(3, dtype='float64'))
        # If we're using a subvolume, we now re-divide.
        if subvolume is not None:
            self._data_source = pf.h.periodic_region_strict([0.] * 3, ds_LE,
                ds_RE)
            # Cut up the volume.
            padded, LE, RE, self._data_source = \
                self.partition_hierarchy_3d(ds=self._data_source,
                padding=0.)
        self.bounds = (LE, RE)
        (LE_padding, RE_padding) = self.padding
        parallelHOPHaloList.__init__(self, self._data_source, self.padding, \
        self.num_neighbors, self.bounds, total_mass, period, \
        threshold=threshold, dm_only=dm_only, rearrange=rearrange,
            premerge=premerge, tree=self.tree)
        self._join_halolists()
        yt_counters("Final Grouping")

    def _subsample_points(self):
        # Read in a random subset of the points in each domain, and then
        # collect them on the root task.
        xp = self._data_source["particle_position_x"]
        n_parts = self.comm.mpi_allreduce(xp.size, op='sum')
        local_parts = xp.size
        random_points = int(self.sample * n_parts)
        # We want to get a representative selection of random particles in
        # each subvolume.
        adjust = float(local_parts) / (float(n_parts) / self.comm.size)
        n_random = int(adjust * float(random_points) / self.comm.size)
        mylog.info("Reading in %d random particles." % n_random)
        # Get unique random particles.
        my_points = na.empty((n_random, 3), dtype='float64')
        uni = na.array(random.sample(xrange(xp.size), n_random))
        uni = uni[uni.argsort()]
        my_points[:, 0] = xp[uni]
        del xp
        self._data_source.clear_data()
        my_points[:, 1] = self._data_source["particle_position_y"][uni]
        self._data_source.clear_data()
        my_points[:, 2] = self._data_source["particle_position_z"][uni]
        self._data_source.clear_data()
        del uni
        # Collect them on the root task.
        mine, sizes = self.comm.mpi_info_dict(n_random)
        if mine == 0:
            tot_random = sum(sizes.values())
            root_points = na.empty((tot_random, 3), dtype='float64')
            root_points.shape = (1, 3 * tot_random)
        else:
            root_points = na.empty([])
        my_points.shape = (1, n_random * 3)
        root_points = self.comm.par_combine_object(my_points[0],
                datatype="array", op="cat")
        del my_points
        if mine == 0:
            root_points.shape = (tot_random, 3)
        return root_points

    def _recursive_divide(self, points, bounds, level, cut_list):
        dim = cut_list[level][0]
        parts = points.shape[0]
        num_bins = 1000
        width = bounds[1][dim] - bounds[0][dim]
        bin_width = width / num_bins
        bins = na.arange(num_bins + 1, dtype='float64') * bin_width + \
            bounds[0][dim]
        counts, bins = na.histogram(points[:, dim], bins)
        # Find the bin that passes the cut points.
        midpoints = [bounds[0][dim]]
        sum = 0
        bin = 0
        for step in xrange(1, cut_list[level][1]):
            while sum < ((parts * step) / cut_list[level][1]):
                lastsum = sum
                sum += counts[bin]
                bin += 1
            # Bin edges
            left_edge = bins[bin - 1]
            right_edge = bins[bin]
            # Find a better approx of the midpoint cut
            # line using a linear approx.
            a = float(sum - lastsum) / (right_edge - left_edge)
            midpoints.append(left_edge + (0.5 - \
                (float(lastsum) / parts / 2)) / a)
        midpoints.append(bounds[1][dim])

        # Split the points & update the bounds.
        subpoints = []
        subbounds = []
        for pair in zip(midpoints[:-1], midpoints[1:]):
            select = na.bitwise_and(points[:, dim] >= pair[0],
                points[:, dim] < pair[1])
            subpoints.append(points[select])
            nb = bounds.copy()
            nb[0][dim] = pair[0]
            nb[1][dim] = pair[1]
            subbounds.append(nb)
        # If we're at the maxlevel, make a bucket. Otherwise, recurse down.
        maxlevel = len(cut_list) - 1
        for pair in zip(subpoints, subbounds):
            if level == maxlevel:
                self.bucket_bounds.append(pair[1])
            else:
                self._recursive_divide(pair[0], pair[1], level + 1, cut_list)

    def _join_halolists(self):
        if self.group_count == 0:
            mylog.info("There are no halos found.")
            return
        ms = -self.Tot_M.copy()
        del self.Tot_M
        Cx = self.CoM[:, 0].copy()
        sorted = na.lexsort([Cx, ms])
        del Cx, ms
        self._groups = self._groups[sorted]
        self._max_dens = self._max_dens[sorted]
        for i in xrange(self.group_count):
            self._groups[i].id = i
        del sorted, self.group_sizes, self.CoM


class HOPHaloFinder(GenericHaloFinder, HOPHaloList):
    def __init__(self, pf, subvolume=None, threshold=160, dm_only=True,
            padding=0.02, total_mass=None):
        r"""HOP halo finder.

        Halos are built by:
        1. Calculating a density for each particle based on a smoothing kernel.
        2. Recursively linking particles to other particles from lower density
        particles to higher.
        3. Geometrically proximate chains are identified and
        4. merged into final halos following merging rules.

        Lower thresholds generally produce more halos, and the largest halos
        become larger. Also, halos become more filamentary and over-connected.

        Eisenstein and Hut. "HOP: A New Group-Finding Algorithm for N-Body
        Simulations." ApJ (1998) vol. 498 pp. 137-142

        Parameters
        ----------
        pf : `StaticOutput`
            The parameter file on which halo finding will be conducted.
        subvolume : `yt.data_objects.api.AMRData`, optional
            A region over which HOP will be run, which can be used to run HOP
            on a subvolume of the full volume. Default = None, which defaults
            to the full volume automatically.
        threshold : float
            The density threshold used when building halos. Default = 160.0.
        dm_only : bool
            If True, only dark matter particles are used when building halos.
            Default = False.
        padding : float
            When run in parallel, the finder needs to surround each subvolume
            with duplicated particles for halo finidng to work. This number
            must be no smaller than the radius of the largest halo in the box
            in code units. Default = 0.02.
        total_mass : float
            If HOP is run on the same dataset mulitple times, the total mass
            of particles in Msun units in the full volume can be supplied here
            to save time.
            This must correspond to the particles being operated on, meaning
            if stars are included in the halo finding, they must be included
            in this mass as well, and visa-versa.
            If halo finding on a subvolume, this still corresponds with the
            mass in the entire volume.
            Default = None, which means the total mass is automatically
            calculated.

        Examples
        --------
        >>> pf = load("RedshiftOutput0000")
        >>> halos = HaloFinder(pf)
        """
        if subvolume is not None:
            ds_LE = na.array(subvolume.left_edge)
            ds_RE = na.array(subvolume.right_edge)
        self.period = pf.domain_right_edge - pf.domain_left_edge
        self._data_source = pf.h.all_data()
        GenericHaloFinder.__init__(self, pf, self._data_source, dm_only,
            padding)
        # do it once with no padding so the total_mass is correct
        # (no duplicated particles), and on the entire volume, even if only
        # a small part is actually going to be used.
        self.padding = 0.0
        padded, LE, RE, self._data_source = \
            self.partition_hierarchy_3d(ds=self._data_source,
                padding=self.padding)
        # For scaling the threshold, note that it's a passthrough
        if total_mass is None:
            if dm_only:
                select = self._get_dm_indices()
                total_mass = \
                    self.comm.mpi_allreduce((self._data_source["ParticleMassMsun"][select]).sum(dtype='float64'), op='sum')
            else:
                total_mass = self.comm.mpi_allreduce(self._data_source.quantities["TotalQuantity"]("ParticleMassMsun")[0], op='sum')
        # MJT: Note that instead of this, if we are assuming that the particles
        # are all on different processors, we should instead construct an
        # object representing the entire domain and sum it "lazily" with
        # Derived Quantities.
        if subvolume is not None:
            self._data_source = pf.h.periodic_region_strict([0.] * 3, ds_LE, ds_RE)
        else:
            self._data_source = pf.h.all_data()
        self.padding = padding  # * pf["unitary"] # This should be clevererer
        padded, LE, RE, self._data_source = \
            self.partition_hierarchy_3d(ds=self._data_source,
            padding=self.padding)
        self.bounds = (LE, RE)
        # sub_mass can be skipped if subvolume is not used and this is not
        # parallel.
        if subvolume is None and \
                ytcfg.getint("yt", "__topcomm_parallel_size") == 1:
            sub_mass = total_mass
        elif dm_only:
            select = self._get_dm_indices()
            sub_mass = self._data_source["ParticleMassMsun"][select].sum(dtype='float64')
        else:
            sub_mass = \
                self._data_source.quantities["TotalQuantity"]("ParticleMassMsun")[0]
        HOPHaloList.__init__(self, self._data_source,
            threshold * total_mass / sub_mass, dm_only)
        self._parse_halolist(total_mass / sub_mass)
        self._join_halolists()


class FOFHaloFinder(GenericHaloFinder, FOFHaloList):
    def __init__(self, pf, subvolume=None, link=0.2, dm_only=True,
        padding=0.02):
        r"""Friends-of-friends halo finder.

        Halos are found by linking together all pairs of particles closer than
        some distance from each other. Particles may have multiple links,
        and halos are found by recursively linking together all such pairs.

        Larger linking lengths produce more halos, and the largest halos
        become larger. Also, halos become more filamentary and over-connected.

        Davis et al. "The evolution of large-scale structure in a universe
        dominated by cold dark matter." ApJ (1985) vol. 292 pp. 371-394

        Parameters
        ----------
        pf : `StaticOutput`
            The parameter file on which halo finding will be conducted.
        subvolume : `yt.data_objects.api.AMRData`, optional
            A region over which HOP will be run, which can be used to run HOP
            on a subvolume of the full volume. Default = None, which defaults
            to the full volume automatically.
        link : float
            If positive, the interparticle distance (compared to the overall
            average) used to build the halos. If negative, this is taken to be
            the *actual* linking length, and no other calculations will be
            applied.  Default = 0.2.
        dm_only : bool
            If True, only dark matter particles are used when building halos.
            Default = False.
        padding : float
            When run in parallel, the finder needs to surround each subvolume
            with duplicated particles for halo finidng to work. This number
            must be no smaller than the radius of the largest halo in the box
            in code units. Default = 0.02.

        Examples
        --------
        >>> pf = load("RedshiftOutput0000")
        >>> halos = FOFHaloFinder(pf)
        """
        if subvolume is not None:
            ds_LE = na.array(subvolume.left_edge)
            ds_RE = na.array(subvolume.right_edge)
        self.period = pf.domain_right_edge - pf.domain_left_edge
        self.pf = pf
        self.hierarchy = pf.h
        self._data_source = pf.h.all_data()
        GenericHaloFinder.__init__(self, pf, self._data_source, dm_only,
            padding)
        self.padding = 0.0  # * pf["unitary"] # This should be clevererer
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self._data_source = \
            self.partition_hierarchy_3d(ds=self._data_source,
            padding=self.padding)
        if link > 0.0:
            n_parts = self.comm.mpi_allreduce(self._data_source["particle_position_x"].size, op='sum')
            # get the average spacing between particles
            #l = pf.domain_right_edge - pf.domain_left_edge
            #vol = l[0] * l[1] * l[2]
            # Because we are now allowing for datasets with non 1-periodicity,
            # but symmetric, vol is always 1.
            vol = 1.
            avg_spacing = (float(vol) / n_parts) ** (1. / 3.)
            linking_length = link * avg_spacing
        else:
            linking_length = na.abs(link)
        self.padding = padding
        if subvolume is not None:
            self._data_source = pf.h.periodic_region_strict([0.] * 3, ds_LE,
                ds_RE)
        else:
            self._data_source = pf.h.all_data()
        padded, LE, RE, self._data_source = \
            self.partition_hierarchy_3d(ds=self._data_source,
            padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        #self._reposition_particles((LE, RE))
        # here is where the FOF halo finder is run
        mylog.info("Using a linking length of %0.3e", linking_length)
        FOFHaloList.__init__(self, self._data_source, linking_length, dm_only)
        self._parse_halolist(1.)
        self._join_halolists()

HaloFinder = HOPHaloFinder


class LoadHaloes(GenericHaloFinder, LoadedHaloList):
    def __init__(self, pf, basename):
        r"""Load the full halo data into memory.

        This function takes the output of `GenericHaloFinder.dump` and
        re-establishes the list of halos in memory. This enables the full set
        of halo analysis features without running the halo finder again. To
        be precise, the particle data for each halo is only read in when
        necessary, so examining a single halo will not require as much memory
        as is required for halo finding.

        Parameters
        ----------
        basename : String
            The base name of the files that will be read in. This should match
            what was used when `GenericHaloFinder.dump` was called. Default =
            "HopAnalysis".

        Examples
        --------
        >>> pf = load("data0005")
        >>> halos = LoadHaloes(pf, "HopAnalysis")
        """
        self.basename = basename
        LoadedHaloList.__init__(self, pf, self.basename)

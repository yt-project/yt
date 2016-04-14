"""
HOP-output data handling



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import gc
from yt.utilities.on_demand_imports import _h5py as h5py
import math
import numpy as np
import os.path as path
from functools import cmp_to_key
from yt.extern.six import add_metaclass
from yt.extern.six.moves import zip as izip

from yt.config import ytcfg
from yt.funcs import mylog, ensure_dir_exists
from yt.utilities.math_utils import \
    get_rotation_matrix, \
    periodic_dist
from yt.utilities.physical_constants import \
    mass_sun_cgs
from yt.utilities.physical_ratios import \
    rho_crit_g_cm3_h2, \
    TINY

from .hop.EnzoHop import RunHOP
from .fof.EnzoFOF import RunFOF

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelDummy, \
    ParallelAnalysisInterface, \
    parallel_blocking_call


@add_metaclass(ParallelDummy)
class Halo(object):
    """
    A data source that returns particle information about the members of a
    HOP-identified halo.
    """
    _distributed = False
    _processing = False
    _owner = 0
    indices = None
    dont_wrap = ["get_sphere", "write_particle_list"]
    extra_wrap = ["__getitem__"]

    def __init__(self, halo_list, id, indices=None, size=None, CoM=None,
        max_dens_point=None, group_total_mass=None, max_radius=None,
        bulk_vel=None, tasks=None, rms_vel=None, supp=None, ptype=None):
        if ptype is None:
            ptype = "all"
        self.ptype = ptype
        self.halo_list = halo_list
        self._max_dens = halo_list._max_dens
        self.id = id
        self.data = halo_list._data_source
        self.ds = self.data.ds
        self.gridsize = (self.ds.domain_right_edge - \
                 self.ds.domain_left_edge)
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
        # A supplementary data dict.
        if supp is None:
            self.supp = {}
        else:
            self.supp = supp
        self._saved_fields = {}
        self._ds_sort = None
        self._particle_mask = None

    @property
    def particle_mask(self):
        # Dynamically create the masking array for particles, and get
        # the data using standard yt methods.
        if self._particle_mask is not None:
            return self._particle_mask
        # This is from disk.
        pid = self.__getitem__('particle_index')
        # This is from the sphere.
        if self._name == "RockstarHalo":
            ds = self.ds.sphere(self.CoM, self._radjust * self.max_radius)
        elif self._name == "LoadedHalo":
            ds = self.ds.sphere(self.CoM, np.maximum(self._radjust * \
            self.ds.quan(self.max_radius, 'code_length'), \
            self.ds.index.get_smallest_dx()))
        sp_pid = ds['particle_index']
        self._ds_sort = sp_pid.argsort()
        sp_pid = sp_pid[self._ds_sort]
        # This matches them up.
        self._particle_mask = np.in1d(sp_pid, pid)
        return self._particle_mask

    def center_of_mass(self):
        r"""Calculate and return the center of mass.

        The center of mass of the halo is directly calculated and returned.

        Examples
        --------
        >>> com = halos[0].center_of_mass()
        """
        if self.CoM is not None:
            return self.CoM
        pm = self["particle_mass"].in_units('Msun')
        c = {}
        # We shift into a box where the origin is the left edge
        c[0] = self["particle_position_x"] - self.ds.domain_left_edge[0]
        c[1] = self["particle_position_y"] - self.ds.domain_left_edge[1]
        c[2] = self["particle_position_z"] - self.ds.domain_left_edge[2]
        com = []
        for i in range(3):
            # A halo is likely periodic around a boundary if the distance
            # between the max and min particle
            # positions are larger than half the box.
            # So skip the rest if the converse is true.
            # Note we might make a change here when periodicity-handling is
            # fully implemented.
            if (c[i].max() - c[i].min()) < (self.ds.domain_width[i] / 2.):
                com.append(c[i])
                continue
            # Now we want to flip around only those close to the left boundary.
            sel = (c[i] <= (self.ds.domain_width[i]/2))
            c[i][sel] += self.ds.domain_width[i]
            com.append(c[i])

        c = (com * pm).sum(axis=1) / pm.sum()
        c = self.ds.arr(c, 'code_length')

        return c%self.ds.domain_width + self.ds.domain_left_edge

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
        return np.array([
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
        return self["particle_mass"].in_units('Msun').sum()

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
        pm = self["particle_mass"].in_units('Msun')
        vx = (self["particle_velocity_x"] * pm).sum()
        vy = (self["particle_velocity_y"] * pm).sum()
        vz = (self["particle_velocity_z"] * pm).sum()
        return self.ds.arr([vx, vy, vz], vx.units) / pm.sum()

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
        pm = self["particle_mass"].in_units('Msun')
        sm = pm.sum()
        vx = (self["particle_velocity_x"] - bv[0]) * pm / sm
        vy = (self["particle_velocity_y"] - bv[1]) * pm / sm
        vz = (self["particle_velocity_z"] - bv[2]) * pm / sm
        s = vx ** 2. + vy ** 2. + vz ** 2.
        ms = np.mean(s)
        return np.sqrt(ms) * pm.size

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
        rx = np.abs(self["particle_position_x"] - center[0])
        ry = np.abs(self["particle_position_y"] - center[1])
        rz = np.abs(self["particle_position_z"] - center[2])
        DW = self.data.ds.domain_right_edge - self.data.ds.domain_left_edge
        r = np.sqrt(np.minimum(rx, DW[0] - rx) ** 2.0
                + np.minimum(ry, DW[1] - ry) ** 2.0
                + np.minimum(rz, DW[2] - rz) ** 2.0)
        return r.max()

    def __getitem__(self, key):
        return self.data[(self.ptype, key)][self.indices]

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
        sphere : `yt.data_objects.api.YTSphere`
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
        sphere = self.data.ds.sphere(center, radius=radius)
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
                   + ["particle_index"]:
            handle.create_dataset("/%s/%s" % (gn, field), data=self[field])
        handle.create_dataset("/%s/particle_mass" % gn,
                              data=self["particle_mass"].in_units('Msun'))
        if ('io','creation_time') in self.data.ds.field_list:
            handle.create_dataset("/%s/creation_time" % gn,
                data=self['creation_time'])
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
        if over.any():
            vir_bin = max(np.arange(bins + 1)[over])
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
        h = self.ds.hubble_constant
        Om_matter = self.ds.omega_matter
        z = self.ds.current_redshift
        period = self.ds.domain_right_edge - \
            self.ds.domain_left_edge
        thissize = self.get_size()
        rho_crit = rho_crit_g_cm3_h2 * h ** 2.0 * Om_matter  # g cm^-3
        Msun2g = mass_sun_cgs
        rho_crit = rho_crit * ((1.0 + z) ** 3.0)
        # Get some pertinent information about the halo.
        self.mass_bins = self.ds.arr(np.zeros(self.bin_count + 1,
                                              dtype='float64'),'Msun')
        dist = np.empty(thissize, dtype='float64')
        cen = self.center_of_mass()
        mark = 0
        # Find the distances to the particles. I don't like this much, but I
        # can't see a way to eliminate a loop like this, either here or in
        # yt.math.
        for pos in izip(self["particle_position_x"],
                self["particle_position_y"], self["particle_position_z"]):
            dist[mark] = periodic_dist(cen, pos, period)
            mark += 1
        # Set up the radial bins.
        # Multiply min and max to prevent issues with digitize below.
        self.radial_bins = np.logspace(math.log10(min(dist) * .99 + TINY),
            math.log10(max(dist) * 1.01 + 2 * TINY), num=self.bin_count + 1)
        self.radial_bins = self.ds.arr(self.radial_bins,'code_length')
        # Find out which bin each particle goes into, and add the particle
        # mass to that bin.
        inds = np.digitize(dist, self.radial_bins) - 1
        if self["particle_position_x"].size > 1:
            for index in np.unique(inds):
                self.mass_bins[index] += \
                np.sum(self["particle_mass"][inds == index]).in_units('Msun')
        # Now forward sum the masses in the bins.
        for i in range(self.bin_count):
            self.mass_bins[i + 1] += self.mass_bins[i]
        # Calculate the over densities in the bins.
        self.overdensity = self.mass_bins * Msun2g / \
            (4./3. * math.pi * rho_crit * \
            (self.radial_bins )**3.0)

    def _get_ellipsoid_parameters_basic(self):
        np.seterr(all='ignore')
        # check if there are 4 particles to form an ellipsoid
        # neglecting to check if the 4 particles in the same plane,
        # that is almost certainly never to occur,
        # will deal with it later if it ever comes up
        if np.size(self["particle_position_x"]) < 4:
            mylog.warning("Too few particles for ellipsoid parameters.")
            return (0, 0, 0, 0, 0, 0, 0)
        # Calculate the parameters that describe the ellipsoid of
        # the particles that constitute the halo. This function returns
        # all the parameters except for the center of mass.
        com = self.center_of_mass()
        position = [self["particle_position_x"],
                    self["particle_position_y"],
                    self["particle_position_z"]]
        # Locate the furthest particle from com, its vector length and index
        DW = np.array([self.gridsize[0],self.gridsize[1],self.gridsize[2]])
        position = [position[0] - com[0],
                    position[1] - com[1],
                    position[2] - com[2]]
        # different cases of particles being on other side of boundary
        for axis in range(np.size(DW)):
            cases = np.array([position[axis],
                                position[axis] + DW[axis],
                              position[axis] - DW[axis]])
            # pick out the smallest absolute distance from com
            position[axis] = np.choose(np.abs(cases).argmin(axis=0), cases)
        # find the furthest particle's index
        r = np.sqrt(position[0]**2 +
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
        e0_vector_copy = np.empty((np.size(position[0]), 3), dtype='float64')
        for i in range(3):
            e0_vector_copy[:, i] = e0_vector[i]
        rr = np.array([position[0],
                       position[1],
                       position[2]]).T # Similar to tB_vector in old code.
        tC_vector = np.cross(e0_vector_copy, rr)
        te2 = tC_vector.copy()
        for dim in range(3):
            te2[:,dim] *= np.sum(tC_vector**2., axis = 1)**(-0.5)
        te1 = np.cross(te2, e0_vector_copy)
        length = np.abs(-np.sum(rr * te1, axis = 1) * \
            (1. - np.sum(rr * e0_vector_copy, axis = 1)**2. * \
            mag_A**-2.)**(-0.5))
        # This problem apparently happens sometimes, that the NaNs are turned
        # into infs, which messes up the nanargmax below.
        length[length == np.inf] = 0.
        tB_index = np.nanargmax(length) # ignores NaNs created above.
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
        length = np.abs(np.sum(rr * temp_e2, axis = 1) * (1 - \
            np.sum(rr * temp_e0, axis = 1)**2. * mag_A**-2. - \
            np.sum(rr * temp_e1, axis = 1)**2. * mag_B**-2.)**(-0.5))
        length[length == np.inf] = 0.
        tC_index = np.nanargmax(length)
        mag_C = length[tC_index]
        # tilt is calculated from the rotation about x axis
        # needed to align e1 vector with the y axis
        # after e0 is aligned with x axis
        # find the t1 angle needed to rotate about z axis to align e0 onto x-z plane
        t1 = np.arctan(-e0_vector[1] / e0_vector[0])
        RZ = get_rotation_matrix(t1, (0, 0, 1))
        r1 = np.dot(RZ, e0_vector)
        # find the t2 angle needed to rotate about y axis to align e0 to x
        t2 = np.arctan(r1[2] / r1[0])
        RY = get_rotation_matrix(t2, (0, 1, 0))
        r2 = np.dot(RY, np.dot(RZ, e1_vector))
        # find the tilt angle needed to rotate about x axis to align e1 to y and e2 to z
        tilt = np.arctan(-r2[2] / r2[1])
        return (mag_A, mag_B, mag_C, e0_vector[0], e0_vector[1],
            e0_vector[2], tilt)

class HOPHalo(Halo):
    _name = "HOPHalo"
    pass


class FOFHalo(Halo):

    def maximum_density(self):
        r"""Not implemented."""
        return -1

    def maximum_density_location(self):
        r"""Not implemented."""
        return self.center_of_mass()


class LoadedHalo(Halo):
    _name = "LoadedHalo"
    # See particle_mask
    _radjust = 1.05

    def __init__(self, ds, id, size=None, CoM=None,
        max_dens_point=None, group_total_mass=None, max_radius=None, bulk_vel=None,
        rms_vel=None, fnames=None, mag_A=None, mag_B=None, mag_C=None,
        e0_vec=None, tilt=None, supp=None):

        self.ds = ds
        self.gridsize = (self.ds.domain_right_edge - \
            self.ds.domain_left_edge)
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
        self.e0_vec = e0_vec
        self.tilt = tilt
        # locs=the names of the h5 files that have particle data for this halo
        self.fnames = fnames
        self.bin_count = None
        self.overdensity = None
        self.indices = np.array([])  # Never used for a LoadedHalo.
        self._saved_fields = {}
        self._ds_sort = None
        self._particle_mask = None
        # A supplementary data dict.
        if supp is None:
            self.supp = {}
        else:
            self.supp = supp
        self._saved_fields = {}
        self._ds_sort = None
        self._particle_mask = None
        self._pid_sort = None


    def __getitem__(self, key):
        # This function will try to get particle data in one of three ways,
        # in descending preference.
        # 1. From saved_fields, e.g. we've already got it.
        # 2. From the halo h5 files off disk.
        # 3. Use the unique particle indexes of the halo to select a missing
        # field from a Sphere.
        if key in self._saved_fields:
            # We've already got it.
            return self._saved_fields[key]
        # Gotta go get it from the halo h5 files.
        field_data = self._get_particle_data(self.id, self.fnames,
            self.size, key)
        if field_data is not None:
            if key == 'particle_index':
                #this is an index for turning data sorted by particle index
                #into the same order as the fields on disk
                self._pid_sort = field_data.argsort().argsort()
            #convert to YTArray using the data from disk
            if key == 'particle_mass':
                field_data = self.ds.arr(field_data, 'Msun')
            else:
                field_data = self.ds.arr(field_data,
                    self.ds._get_field_info('unknown',key).units)
            self._saved_fields[key] = field_data
            return self._saved_fields[key]
        # We won't store this field below in saved_fields because
        # that would mean keeping two copies of it, one in the yt
        # machinery and one here.
        ds = self.ds.sphere(self.CoM, np.maximum(self._radjust * \
            self.ds.quan(self.max_radius, 'code_length'), \
            self.ds.index.get_smallest_dx()))
        # If particle_mask hasn't been called once then _ds_sort won't have
        # the proper values set yet
        if self._particle_mask is None:
            self.particle_mask
        return ds[key][self._ds_sort][self.particle_mask][self._pid_sort]

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
            field_data = np.empty(size, dtype='int64')
        else:
            field_data = np.empty(size, dtype='float64')
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
            return (self.mag_A, self.mag_B, self.mag_C, self.e0_vec[0],
                self.e0_vec[1], self.e0_vec[2], self.tilt)
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

        basic_parameters = self._get_ellipsoid_parameters_basic_loadedhalo()
        toreturn = [self.center_of_mass()]
        updated = [basic_parameters[0], basic_parameters[1],
            basic_parameters[2], np.array([basic_parameters[3],
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
        ellipsoid : `yt.data_objects.data_containers.YTEllipsoid`
            The ellipsoidal data object.

        Examples
        --------
        >>> ell = halos[0].get_ellipsoid()
        """
        ep = self.get_ellipsoid_parameters()
        ell = self.ds.ellipsoid(ep[0], ep[1], ep[2], ep[3], ep[4], ep[5])
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
        sphere : `yt.data_objects.api.YTSphere`
            The empty data source.

        Examples
        --------
        >>> sp = halos[0].get_sphere()
        """
        cen = self.center_of_mass()
        r = self.maximum_radius()
        return self.ds.sphere(cen, r)

class TextHalo(LoadedHalo):
    def __init__(self, ds, id, size=None, CoM=None,

        max_dens_point=None, group_total_mass=None, max_radius=None, bulk_vel=None,
        rms_vel=None, fnames=None, mag_A=None, mag_B=None, mag_C=None,
        e0_vec=None, tilt=None, supp=None):

        self.ds = ds
        self.gridsize = (self.ds.domain_right_edge - \
            self.ds.domain_left_edge)
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
        self.e0_vec = e0_vec
        self.tilt = tilt
        self.bin_count = None
        self.overdensity = None
        self.indices = np.array([])  # Never used for a LoadedHalo.
        # A supplementary data dict.
        if supp is None:
            self.supp = {}
        else:
            self.supp = supp

    def __getitem__(self, key):
        # We'll just pull it from the sphere.
        return self.get_sphere()[key]

    def maximum_density(self):
        r"""Undefined for text halos."""
        return -1

    def maximum_density_location(self):
        r"""Undefined, default to CoM"""
        return self.center_of_mass()

    def get_size(self):
        # Have to just get it from the sphere.
        return self["particle_position_x"].size


class HaloList(object):

    _fields = ["particle_position_%s" % ax for ax in 'xyz']

    def __init__(self, data_source, dm_only=True, redshift=-1,
                 ptype=None):
        """
        Run hop on *data_source* with a given density *threshold*.  If
        *dm_only* is True (default), only run it on the dark matter particles,
        otherwise on all particles.  Returns an iterable collection of
        *HopGroup* items.
        """
        self._data_source = data_source
        self.dm_only = dm_only
        if ptype is None:
            ptype = "all"
        self.ptype = ptype
        self._groups = []
        self._max_dens = {}
        self.__obtain_particles()
        self._run_finder()
        mylog.info("Parsing outputs")
        self._parse_output()
        mylog.debug("Finished. (%s)", len(self))
        self.redshift = redshift

    def __obtain_particles(self):
        if self.dm_only:
            ii = self._get_dm_indices()
        else:
            ii = slice(None)
        self.particle_fields = {}
        for field in self._fields:
            tot_part = self._data_source[(self.ptype, field)].size
            if field == "particle_index":
                self.particle_fields[field] = \
                    self._data_source[(self.ptype, field)][ii].astype('int64')
            else:
                self.particle_fields[field] = \
                    self._data_source[(self.ptype, field)][ii].astype('float64')
            del self._data_source[(self.ptype, field)]
        self._base_indices = np.arange(tot_part)[ii]
        gc.collect()

    def _get_dm_indices(self):
        if ('io','creation_time') in self._data_source.index.field_list:
            mylog.debug("Differentiating based on creation time")
            return (self._data_source["creation_time"] <= 0)
        elif ('io','particle_type') in self._data_source.index.field_list:
            mylog.debug("Differentiating based on particle type")
            return (self._data_source["particle_type"] == 1)
        else:
            mylog.warning("No particle_type, no creation_time, so not distinguishing.")
            return slice(None)

    def _parse_output(self):
        unique_ids = np.unique(self.tags)
        counts = np.bincount(self.tags + 1)
        sort_indices = np.argsort(self.tags)
        grab_indices = np.indices(self.tags.shape).ravel()[sort_indices]
        dens = self.densities[sort_indices]
        cp = 0
        for i in unique_ids:
            cp_c = cp + counts[i + 1]
            if i == -1:
                cp += counts[i + 1]
                continue
            group_indices = grab_indices[cp:cp_c]
            self._groups.append(self._halo_class(self, i, group_indices,
                                                 ptype=self.ptype))
            md_i = np.argmax(dens[cp:cp_c])
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
        f.write("# REDSHIFT OF OUTPUT = %f\n" % (self.redshift))

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
                               "mag_A", "mag_B", "mag_C", "e0_vec0",
                               "e0_vec1", "e0_vec2", "tilt", "\n"]))

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
    """
    Run hop on *data_source* with a given density *threshold*.  If
    *dm_only* is True (default), only run it on the dark matter particles, otherwise
    on all particles.  Returns an iterable collection of *HopGroup* items.
    """
    _name = "HOP"
    _halo_class = HOPHalo
    _fields = ["particle_position_%s" % ax for ax in 'xyz'] + \
              ["particle_mass"]

    def __init__(self, data_source, threshold=160.0, dm_only=True,
                 ptype=None):
        self.threshold = threshold
        mylog.info("Initializing HOP")
        HaloList.__init__(self, data_source, dm_only, ptype=ptype)

    def _run_finder(self):
        self.densities, self.tags = \
            RunHOP(self.particle_fields["particle_position_x"] / self.period[0],
                self.particle_fields["particle_position_y"] / self.period[1],
                self.particle_fields["particle_position_z"] / self.period[2],
                self.particle_fields["particle_mass"].in_units('Msun'),
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

    def __init__(self, data_source, link=0.2, dm_only=True, redshift=-1,
                 ptype=None):
        self.link = link
        mylog.info("Initializing FOF")
        HaloList.__init__(self, data_source, dm_only, redshift=redshift,
                          ptype=ptype)

    def _run_finder(self):
        self.tags = \
        RunFOF(self.particle_fields["particle_position_x"] / self.period[0],
               self.particle_fields["particle_position_y"] / self.period[1],
               self.particle_fields["particle_position_z"] / self.period[2],
               self.link)
        self.densities = np.ones(self.tags.size, dtype='float64') * -1
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

    def __init__(self, ds, basename):
        ParallelAnalysisInterface.__init__(self)
        self.ds = ds
        self._groups = []
        self.basename = basename
        self._retrieve_halos()

    def _retrieve_halos(self):
        # First get the halo particulars.
        with open("%s.out" % self.basename, 'r') as fh:
            lines = fh.readlines()
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
            CoM = np.array([float(line[7]), float(line[8]), float(line[9])])
            max_dens_point = np.array([float(line[3]), float(line[4]),
                float(line[5]), float(line[6])])
            group_total_mass = float(line[1])
            max_radius = float(line[13])
            bulk_vel = np.array([float(line[10]), float(line[11]),
                float(line[12])])
            rms_vel = float(line[14])
            if len(line) == 15:
                # No ellipsoid information
                self._groups.append(LoadedHalo(self.ds, halo, size = size,
                    CoM = CoM,
                    max_dens_point = max_dens_point,
                    group_total_mass = group_total_mass, max_radius = max_radius,
                    bulk_vel = bulk_vel, rms_vel = rms_vel, fnames = fnames))
            elif len(line) == 22:
                # Ellipsoid information
                mag_A = float(line[15])
                mag_B = float(line[16])
                mag_C = float(line[17])
                e0_vec0 = float(line[18])
                e0_vec1 = float(line[19])
                e0_vec2 = float(line[20])
                e0_vec = np.array([e0_vec0, e0_vec1, e0_vec2])
                tilt = float(line[21])
                self._groups.append(LoadedHalo(self.ds, halo, size = size,
                    CoM = CoM,
                    max_dens_point = max_dens_point,
                    group_total_mass = group_total_mass, max_radius = max_radius,
                    bulk_vel = bulk_vel, rms_vel = rms_vel, fnames = fnames,
                    mag_A = mag_A, mag_B = mag_B, mag_C = mag_C, e0_vec = e0_vec,
                    tilt = tilt))
            else:
                mylog.error("I am unable to parse this line. Too many or too few items. %s" % orig)
            halo += 1

    def _collect_halo_data_locations(self):
        # The halos are listed in order in the file.
        with open("%s.txt" % self.basename, 'r') as fh:
            lines = fh.readlines()
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
        return locations

class TextHaloList(HaloList):
    _name = "Text"

    def __init__(self, ds, fname, columns, comment):
        ParallelAnalysisInterface.__init__(self)
        self.ds = ds
        self._groups = []
        self._retrieve_halos(fname, columns, comment)

    def _retrieve_halos(self, fname, columns, comment):
        # First get the halo particulars.
        with open(fname, 'r') as fh:
            lines = fh.readlines()
        halo = 0
        base_set = ['x', 'y', 'z', 'r']
        keys = columns.keys()
        extra = (len(keys) > 4)
        for line in lines:
            # Skip commented lines.
            if line[0] == comment: continue
            line = line.split()
            x = float(line[columns['x']])
            y = float(line[columns['y']])
            z = float(line[columns['z']])
            r = float(line[columns['r']])
            cen = np.array([x, y, z])
            # Now we see if there's anything else.
            if extra:
                temp_dict = {}
                for key in columns:
                    if key not in base_set:
                        val = float(line[columns[key]])
                        temp_dict[key] = val
            self._groups.append(TextHalo(self.ds, halo,
                CoM = cen, max_radius = r, supp = temp_dict))
            halo += 1

class GenericHaloFinder(HaloList, ParallelAnalysisInterface):
    def __init__(self, ds, data_source, padding=0.0, ptype=None):
        ParallelAnalysisInterface.__init__(self)
        self.ds = ds
        self.index = ds.index
        self.center = (np.array(data_source.right_edge) +
                       np.array(data_source.left_edge)) / 2.0
        if ptype is None:
            ptype = "all"
        self.ptype = ptype

    def _parse_halolist(self, threshold_adjustment):
        groups = []
        max_dens = {}
        hi = 0
        LE, RE = self.bounds
        for halo in self._groups:
            this_max_dens = halo.maximum_density_location()
            # if the most dense particle is in the box, keep it
            if np.all((this_max_dens >= LE) & (this_max_dens <= RE)):
                # Now we add the halo information to OURSELVES, taken from the
                # self.hop_list
                # We need to mock up the HOPHaloList thingie, so we need to
                #     set self._max_dens
                max_dens_temp = list(self._max_dens[halo.id])[0] / \
                    threshold_adjustment
                max_dens[hi] = [max_dens_temp] + \
                    list(self._max_dens[halo.id])[1:4]
                groups.append(self._halo_class(self, hi, ptype=self.ptype))
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
        self._groups = [self._halo_class(self, i, ptype=self.ptype)
                        for i in range(my_first_id)] + \
                       self._groups + \
                       [self._halo_class(self, i, ptype=self.ptype)
                        for i in range(after, nhalos)]
        id = 0
        for proc in sorted(halo_info.keys()):
            for halo in self._groups[id:id + halo_info[proc]]:
                halo.id = id
                halo._distributed = self._distributed
                halo._owner = proc
                id += 1

        def haloCmp(h1, h2):
            def cmp(a, b):
                return (a > b) - (a < b)
            c = cmp(h1.total_mass(), h2.total_mass())
            if c != 0:
                return -1 * c
            if c == 0:
                return cmp(h1.center_of_mass()[0], h2.center_of_mass()[0])
        self._groups.sort(key=cmp_to_key(haloCmp))
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
        dw = self.ds.domain_right_edge - self.ds.domain_left_edge
        for i, ax in enumerate('xyz'):
            arr = self._data_source[self.ptype, "particle_position_%s" % ax]
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
        ensure_dir_exists(filename)
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
        ensure_dir_exists(prefix)
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
        ensure_dir_exists(prefix)
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
        ensure_dir_exists(basename)
        self.write_out("%s.out" % basename, ellipsoid_data)
        self.write_particle_lists(basename)
        self.write_particle_lists_txt(basename)

class HOPHaloFinder(GenericHaloFinder, HOPHaloList):
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
    ds : `Dataset`
        The dataset on which halo finding will be conducted.
    subvolume : `yt.data_objects.data_containers.YTSelectionContainer`, optional
        A region over which HOP will be run, which can be used to run HOP
        on a subvolume of the full volume. Default = None, which defaults
        to the full volume automatically.
    threshold : float
        The density threshold used when building halos. Default = 160.0.
    dm_only : bool (deprecated)
        If True, only dark matter particles are used when building halos.
        This has been deprecated.  Instead, the ptype keyword should be
        used to specify a particle type.
        Default = True.
    ptype : string
        When dm_only is set to False, this sets the type of particle to be
        used for halo finding, with a default of "all".  This should not be
        used when dm_only is set to True.
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
    >>> ds = load("RedshiftOutput0000")
    >>> halos = HaloFinder(ds)
    """
    def __init__(self, ds, subvolume=None, threshold=160, dm_only=True,
                 ptype=None, padding=0.02, total_mass=None):
        if subvolume is not None:
            ds_LE = np.array(subvolume.left_edge)
            ds_RE = np.array(subvolume.right_edge)
        self.period = ds.domain_right_edge - ds.domain_left_edge
        self._data_source = ds.all_data()
        GenericHaloFinder.__init__(self, ds, self._data_source, padding,
                                   ptype=ptype)
        # do it once with no padding so the total_mass is correct
        # (no duplicated particles), and on the entire volume, even if only
        # a small part is actually going to be used.
        self.padding = 0.0
        padded, LE, RE, self._data_source = \
            self.partition_index_3d(ds=self._data_source,
                padding=self.padding)

        if dm_only:
            mylog.warn("dm_only is deprecated.  " +
                       "Use ptype to specify a particle type, instead.")

        # Don't allow dm_only=True and setting a ptype.
        if dm_only and ptype is not None:
            raise RuntimeError(
                "If dm_only is True, ptype must be None.  " + \
                "dm_only must be False if ptype is set.")

        if ptype is None:
            ptype = "all"
        self.ptype = ptype

        # For scaling the threshold, note that it's a passthrough
        if total_mass is None:
            if dm_only:
                select = self._get_dm_indices()
                total_mass = \
                    self.comm.mpi_allreduce((self._data_source['all', "particle_mass"][select].in_units('Msun')).sum(dtype='float64'), op='sum')
            else:
                total_mass = self.comm.mpi_allreduce(
                    self._data_source.quantities.total_quantity(
                        (self.ptype, "particle_mass")).in_units('Msun'), op='sum')
        # MJT: Note that instead of this, if we are assuming that the particles
        # are all on different processors, we should instead construct an
        # object representing the entire domain and sum it "lazily" with
        # Derived Quantities.
        if subvolume is not None:
            self._data_source = ds.region([0.] * 3, ds_LE, ds_RE)
        else:
            self._data_source = ds.all_data()
        self.padding = padding  # * ds["unitary"] # This should be clevererer
        padded, LE, RE, self._data_source = \
            self.partition_index_3d(ds=self._data_source,
            padding=self.padding)
        self.bounds = (LE, RE)
        # sub_mass can be skipped if subvolume is not used and this is not
        # parallel.
        if subvolume is None and \
                ytcfg.getint("yt", "__topcomm_parallel_size") == 1:
            sub_mass = total_mass
        elif dm_only:
            select = self._get_dm_indices()
            sub_mass = self._data_source["particle_mass"][select].in_units('Msun').sum(dtype='float64')
        else:
            sub_mass = \
                self._data_source.quantities.total_quantity(
                    (self.ptype, "particle_mass")).in_units('Msun')
        HOPHaloList.__init__(self, self._data_source,
            threshold * total_mass / sub_mass, dm_only, ptype=self.ptype)
        self._parse_halolist(total_mass / sub_mass)
        self._join_halolists()


class FOFHaloFinder(GenericHaloFinder, FOFHaloList):
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
    ds : `Dataset`
        The dataset on which halo finding will be conducted.
    subvolume : `yt.data_objects.data_containers.YTSelectionContainer`, optional
        A region over which HOP will be run, which can be used to run HOP
        on a subvolume of the full volume. Default = None, which defaults
        to the full volume automatically.
    link : float
        If positive, the interparticle distance (compared to the overall
        average) used to build the halos. If negative, this is taken to be
        the *actual* linking length, and no other calculations will be
        applied.  Default = 0.2.
    dm_only : bool (deprecated)
        If True, only dark matter particles are used when building halos.
        This has been deprecated.  Instead, the ptype keyword should be
        used to specify a particle type.
        Default = True.
    ptype : string
        When dm_only is set to False, this sets the type of particle to be
        used for halo finding, with a default of "all".  This should not be
        used when dm_only is set to True.
    padding : float
        When run in parallel, the finder needs to surround each subvolume
        with duplicated particles for halo finidng to work. This number
        must be no smaller than the radius of the largest halo in the box
        in code units. Default = 0.02.

    Examples
    --------
    >>> ds = load("RedshiftOutput0000")
    >>> halos = FOFHaloFinder(ds)
    """
    def __init__(self, ds, subvolume=None, link=0.2, dm_only=True,
                 ptype=None, padding=0.02):
        if subvolume is not None:
            ds_LE = np.array(subvolume.left_edge)
            ds_RE = np.array(subvolume.right_edge)
        self.period = ds.domain_right_edge - ds.domain_left_edge
        self.ds = ds
        self.index = ds.index
        self.redshift = ds.current_redshift
        self._data_source = ds.all_data()
        GenericHaloFinder.__init__(self, ds, self._data_source, padding)
        self.padding = 0.0  # * ds["unitary"] # This should be clevererer
        # get the total number of particles across all procs, with no padding
        padded, LE, RE, self._data_source = \
            self.partition_index_3d(ds=self._data_source,
            padding=self.padding)

        if dm_only:
            mylog.warn("dm_only is deprecated.  " +
                       "Use ptype to specify a particle type, instead.")

        # Don't allow dm_only=True and setting a ptype.
        if dm_only and ptype is not None:
            raise RuntimeError(
                "If dm_only is True, ptype must be None.  " + \
                "dm_only must be False if ptype is set.")

        if ptype is None:
            ptype = "all"
        self.ptype = ptype

        if link > 0.0:
            n_parts = self.comm.mpi_allreduce(self._data_source["particle_position_x"].size, op='sum')
            # get the average spacing between particles
            #l = ds.domain_right_edge - ds.domain_left_edge
            #vol = l[0] * l[1] * l[2]
            # Because we are now allowing for datasets with non 1-periodicity,
            # but symmetric, vol is always 1.
            vol = 1.
            avg_spacing = (float(vol) / n_parts) ** (1. / 3.)
            linking_length = link * avg_spacing
        else:
            linking_length = np.abs(link)
        self.padding = padding
        if subvolume is not None:
            self._data_source = ds.region([0.] * 3, ds_LE,
                ds_RE)
        else:
            self._data_source = ds.all_data()
        padded, LE, RE, self._data_source = \
            self.partition_index_3d(ds=self._data_source,
            padding=self.padding)
        self.bounds = (LE, RE)
        # reflect particles around the periodic boundary
        #self._reposition_particles((LE, RE))
        # here is where the FOF halo finder is run
        mylog.info("Using a linking length of %0.3e", linking_length)
        FOFHaloList.__init__(self, self._data_source, linking_length, dm_only,
                             redshift=self.redshift, ptype=self.ptype)
        self._parse_halolist(1.)
        self._join_halolists()

HaloFinder = HOPHaloFinder


class LoadHaloes(GenericHaloFinder, LoadedHaloList):
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
    >>> ds = load("data0005")
    >>> halos = LoadHaloes(ds, "HopAnalysis")
    """
    def __init__(self, ds, basename):
        self.basename = basename
        LoadedHaloList.__init__(self, ds, self.basename)

class LoadTextHaloes(GenericHaloFinder, TextHaloList):
    r"""Load a text file of halos.

    Like LoadHaloes, but when all that is available is a plain
    text file. This assumes the text file has the 3-positions of halos
    along with a radius. The halo objects created are spheres.

    Parameters
    ----------
    fname : String
        The name of the text file to read in.

    columns : dict
        A dict listing the column name : column number pairs for data
        in the text file. It is zero-based (like Python).
        An example is {'x':0, 'y':1, 'z':2, 'r':3, 'm':4}.
        Any column name outside of ['x', 'y', 'z', 'r'] will be attached
        to each halo object in the supplementary dict 'supp'. See
        example.

    comment : String
        If the first character of a line is equal to this, the line is
        skipped. Default = "#".

    Examples
    --------
    >>> ds = load("data0005")
    >>> halos = LoadTextHaloes(ds, "list.txt",
        {'x':0, 'y':1, 'z':2, 'r':3, 'm':4},
        comment = ";")
    >>> halos[0].supp['m']
        3.28392048e14
    """
    def __init__(self, ds, filename, columns, comment = "#"):
        TextHaloList.__init__(self, ds, filename, columns, comment)

LoadTextHalos = LoadTextHaloes

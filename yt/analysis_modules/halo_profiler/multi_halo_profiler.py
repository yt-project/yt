"""
HaloProfiler class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

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

import numpy as na
import os
import h5py
import types

from yt.funcs import *
from yt.utilities.math_utils import periodic_dist

from yt.convenience import \
    load
from yt.data_objects.profiles import \
    BinnedProfile1D, EmptyProfileData
from yt.analysis_modules.halo_finding.api import *
from .halo_filters import \
    VirialFilter
from .centering_methods import \
    centering_registry
from yt.data_objects.field_info_container import \
    add_field

from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, \
    parallel_blocking_call, \
    parallel_root_only, \
    parallel_objects
from yt.visualization.fixed_resolution import \
    FixedResolutionBuffer
from yt.visualization.image_writer import write_image

PROFILE_RADIUS_THRESHOLD = 2

class HaloProfiler(ParallelAnalysisInterface):
    "Radial profiling, filtering, and projections for halos in cosmological simulations."
    def __init__(self, dataset, output_dir=None,
                 halos='multiple', halo_list_file='HopAnalysis.out', 
                 halo_list_format='yt_hop', halo_finder_function=parallelHF, 
                 halo_finder_args=None, 
                 halo_finder_kwargs=dict(threshold=160.0, safety=1.5, 
                                         dm_only=False, resize=True, 
                                         fancy_padding=True, rearrange=True),
                 halo_radius=None, radius_units='1', n_profile_bins=50,
                 recenter = None,
                 profile_output_dir='radial_profiles', projection_output_dir='projections',
                 projection_width=8.0, projection_width_units='mpc', project_at_level='max',
                 velocity_center=['bulk', 'halo'], filter_quantities=['id', 'center', 'r_max'], 
                 use_critical_density=False):
        r"""Initialize a Halo Profiler object.
        
        In order to run the halo profiler, the Halo Profiler object must be
        instantiated. At the minimum, the path to a parameter file
        must be provided as the first term.
        
        Parameters
        ----------
        
        dataset : string, required
            The path to the parameter file for the dataset to be analyzed.
        output_dir : string, optional
            If specified, all output will be put into this path instead of 
            in the dataset directories.  Default: None.
        halos :  {"multiple", "single"}, optional
            For profiling more than one halo.  In this mode halos are read in 
            from a list or identified with a halo finder.  In "single" mode,
            the one and only halo 
            center is identified automatically as the location of the peak
            in the density field.  
            Default: "multiple".
        halo_list_file : string, optional
            The name of a file containing the list of halos.  The HaloProfiler
            will  look for this file in the data directory.
            Default: "HopAnalysis.out".
        halo_list_format : {string, dict}
            The format of the halo list file.  "yt_hop" for the format 
            given by yt's halo finders.  "enzo_hop" for the format written
            by enzo_hop. "p-groupfinder"  for P-Groupfinder.  This keyword 
            can also be given in the form of a dictionary specifying the
            column in which various properties can be found.
            For example, {"id": 0, "center": [1, 2, 3], "mass": 4, "radius": 5}.
            Default: "yt_hop".
        halo_finder_function : function
            If halos is set to multiple and the file given by 
            halo_list_file does not exit, the halo finding function
            specified here will be called.  
            Default: HaloFinder (yt_hop).
        halo_finder_args : tuple
            Args given with call to halo finder function.  Default: None.
        halo_finder_kwargs : dict
            kwargs given with call to halo finder function. Default: None.
        recenter : {string, function
            The name of a function that recenters the halo for analysis.
            Default: None.
        halo_radius : float
            If no halo radii are provided in the halo list file, this
            parameter is used to specify the radius out to which radial
            profiles will be made.  This keyword is also 
            used when halos is set to single.  Default: 0.1.
        radius_units : string
            The units of halo_radius.  Default: "1" (code units).
        n_profile_bins : int
            The number of bins in the radial profiles.  Default: 50.
        profile_output_dir : str
            The subdirectory, inside the data directory, in which radial profile 
            output files will be created.
            The directory will be created if it does not exist.  
            Default: "radial_profiles".
        projection_output_dir : str
            The subdirectory, inside the data directory, in which projection 
            output files will be created.
            The directory will be created if it does not exist.  
            Default: "projections".
        projection_width : float
            The width of halo projections.  Default: 8.0.
        projection_width_units : string
            The units of projection_width. Default: "mpc".
        project_at_level : {"max", int}
            The maximum refinement level to be included in projections.  
            Default: "max" (maximum level within the dataset).
        velocity_center  : array_like
            The method in which the halo bulk velocity is calculated (used for 
            calculation of radial and tangential velocities.  Valid options are:
     	        * ["bulk", "halo"] (Default): the velocity provided in
     	          the halo list
                * ["bulk", "sphere"]: the bulk velocity of the sphere
                  centered on the halo center.
    	        * ["max", field]: the velocity of the cell that is the
    	          location of the maximum of the field 
                  specified (used only when halos set to single).
        filter_quantities : array_like
            Quantities from the original halo list file to be written out in the 
            filtered list file.  Default: ['id','center'].
        use_critical_density : bool
            If True, the definition of overdensity for virial quantities
            is calculated with respect to the critical density.
            If False, overdensity is with respect to mean matter density,
            which is lower by a factor of Omega_M.  Default: False.
        
        Examples
        --------
        >>> import yt.analysis_modules.halo_profiler.api as HP
        >>> hp = HP.halo_profiler("DD0242/DD0242")
        
        """
        ParallelAnalysisInterface.__init__(self)

        self.dataset = dataset
        self.output_dir = output_dir
        self.profile_output_dir = profile_output_dir
        self.projection_output_dir = projection_output_dir
        self.n_profile_bins = n_profile_bins
        self.projection_width = projection_width
        self.projection_width_units = projection_width_units
        self.project_at_level = project_at_level
        self.filter_quantities = filter_quantities
        if self.filter_quantities is None: self.filter_quantities = []
        self.use_critical_density = use_critical_density

        self.profile_fields = []
        self.projection_fields = []

        self._halo_filters = []
        self.all_halos = []
        self.filtered_halos = []

        # Create output directory if specified
        if self.output_dir is not None:
            self.__check_directory(self.output_dir)

        # Set halo finder function and parameters, if needed.
        self.halo_finder_function = halo_finder_function
        self.halo_finder_args = halo_finder_args
        if self.halo_finder_args is None: self.halo_finder_args = ()
        self.halo_finder_kwargs = halo_finder_kwargs
        if self.halo_finder_kwargs is None: self.halo_finder_kwargs = {}

        # Set option to get halos from hop or single halo at density maximum.
        # multiple: get halos from hop
        # single: get single halo from density maximum
        self.halos = halos
        if not(self.halos is 'multiple' or self.halos is 'single'):
            mylog.error("Keyword, halos, must be either 'single' or 'multiple'.")
            return None

        # Set halo list format.
        # 'yt_hop': yt hop output.
        # 'enzo_hop': enzo_hop output.
        # dictionary: a dictionary containing fields and their corresponding columns.
        self.halo_list_file = halo_list_file
        if halo_list_format == 'yt_hop':
            self.halo_list_format = {'id':0, 'mass':1, 'np': 2, 
                                     'center':[7, 8, 9], 'velocity':[10, 11, 12], 'r_max':13}
        elif halo_list_format == 'enzo_hop':
            self.halo_list_format = {'id':0, 'center':[4, 5, 6]}
        elif halo_list_format == 'p-groupfinder':
            self.halo_list_format = {'id':3, 'mass':5, 'center':[0, 1, 2], 'r200kpc':8}
        elif isinstance(halo_list_format, types.DictType):
            self.halo_list_format = halo_list_format
        else:
            mylog.error("Keyword, halo_list_format, must be 'yt_hop', 'enzo_hop', 'p-groupfinder', or a dictionary of custom settings.")
            return None

        # Option to recenter sphere someplace else.
        self.recenter = recenter

        # Look for any field that might need to have the bulk velocity set.
        self._need_bulk_velocity = False
        for field in [hp['field'] for hp in self.profile_fields]:
            if 'Velocity' in field or 'Mach' in field:
                self._need_bulk_velocity = True
                break

        # Check validity for VelocityCenter parameter which toggles how the 
        # velocity is zeroed out for radial velocity profiles.
        self.velocity_center = velocity_center[:]
        if self.velocity_center[0] == 'bulk':
            if self.velocity_center[1] == 'halo' and \
                    self.halos is 'single':
                mylog.error("Parameter, VelocityCenter, must be set to 'bulk sphere' or 'max <field>' with halos flag set to 'single'.")
                return None
            if self.velocity_center[1] == 'halo' and \
                    self.halo_list_format is 'enzo_hop':
                mylog.error("Parameter, VelocityCenter, must be 'bulk sphere' for old style hop output files.")
                return None
            if not(self.velocity_center[1] == 'halo' or 
                   self.velocity_center[1] == 'sphere'):
                mylog.error("Second value of VelocityCenter must be either 'halo' or 'sphere' if first value is 'bulk'.")
                return None
        elif self.velocity_center[0] == 'max':
            if self.halos is 'multiple':
                mylog.error("Getting velocity center from a max field value only works with halos='single'.")
                return None
        else:
            mylog.error("First value of parameter, VelocityCenter, must be either 'bulk' or 'max'.")
            return None

        # Create dataset object.
        self.pf = load(self.dataset)
        self.pf.h

        # Figure out what max radius to use for profiling.
        if halo_radius is not None:
            self.halo_radius = halo_radius / self.pf[radius_units]
        elif self.halos is 'single' or not 'r_max' in self.halo_list_format:
            self.halo_radius = 0.1
        else:
            self.halo_radius = None

        # Get halo(s).
        if self.halos is 'single':
            v, center = self.pf.h.find_max('Density')
            singleHalo = {}
            singleHalo['center'] = center
            singleHalo['r_max'] = self.halo_radius * self.pf.units['mpc']
            singleHalo['id'] = 0
            self.all_halos.append(singleHalo)
        elif self.halos is 'multiple':
            # Get hop data.
            self._load_halo_data()
            if len(self.all_halos) == 0:
                mylog.error("No halos loaded, there will be nothing to do.")
                return None
        else:
            mylog.error("I don't know whether to get halos from hop or from density maximum.  This should not have happened.")
            return None

    def add_halo_filter(self, function, *args, **kwargs):
        r"""Filters can be added to create a refined list of halos based on
        their profiles or to avoid profiling halos altogether based on
        information given in the halo list file.
        
        It is often the case that one is looking to identify halos with a
        specific set of properties. This can be accomplished through the
        creation of filter functions. A filter function can take as many args
        and kwargs as you like, as long as the first argument is a profile
        object, or at least a dictionary which contains the profile arrays
        for each field. Filter functions must return a list of two things.
        The first is a True or False indicating whether the halo passed the
        filter. The second is a dictionary containing quantities calculated 
        for that halo that will be written to a file if the halo passes the
        filter. A sample filter function based on virial quantities can be
        found in yt/analysis_modules/halo_profiler/halo_filters.py.
        
        Parameters
        ----------
        function : function
            The name of a halo filter function.
        args : values
            Arguments passed to the halo filter function.
        kwargs : values
            Arguments passed to the halo filter function.
        
        Examples
        -------
        >>> hp.add_halo_filter(HP.VirialFilter, must_be_virialized=True,
                overdensity_field='ActualOverdensity',
                virial_overdensity=200,
                virial_filters=[['TotalMassMsun','>=','1e14']],
                virial_quantities=['TotalMassMsun','RadiusMpc'])
        
        """

        self._halo_filters.append({'function':function, 'args':args, 'kwargs':kwargs})

    def add_profile(self, field, weight_field=None, accumulation=False):
        r"""Add a field for profiling.
        
        Once the halo profiler object has been instantiated,
        fields can be added for profiling using this function. This function
        may be called multiple times, once per field to be added.
        
        Parameters
        ----------
        field : string
            The name of the field.
        weight_field : {None, string}, optional
            The field that will be used to weight the field `field` when
            the radial binning is done. Default: None.
        accumulation : bool
            Whether or not the `field` values should be summed up with the
            radius of the profile.
        
        Examples
        >>> hp.add_profile('CellVolume', weight_field=None, accumulation=True)
        >>> hp.add_profile('TotalMassMsun', weight_field=None, accumulation=True)
        >>> hp.add_profile('Density', weight_field=None, accumulation=False)
        >>> hp.add_profile('Temperature', weight_field='CellMassMsun', accumulation=False)
            
        """

        self.profile_fields.append({'field':field, 'weight_field':weight_field, 
                                    'accumulation':accumulation})

    def add_projection(self, field, weight_field=None, cmap='algae'):
        r"""Make a projection of the specified field.
        
        For the given field, a projection will be produced that can be saved
        to HDF5 or image format. See `make_projections`.
        
        Parameters
        ----------
        field : string
            The name of the field.
        weight_field : string
            The field that will be used to weight the field `field` when
            the projection is done. Default: None.
        cmap : string
            The name of the matplotlib color map that will be used if an
            image is made from the projection. Default="algae".
        
        Examples
        --------
        >>> hp.add_projection('Density', weight_field=None)
        >>> hp.add_projection('Temperature', weight_field='Density')
        >>> hp.add_projection('Metallicity', weight_field='Density')

        """

        self.projection_fields.append({'field':field, 'weight_field':weight_field, 
                                       'cmap': cmap})

    @parallel_blocking_call
    def make_profiles(self, filename=None, prefilters=None, **kwargs):
        r"""Make radial profiles for all halos in the list.
        
        After all the calls to `add_profile`, this will trigger the actual
        calculations and output the profiles to disk.
        
        Paramters
        ---------
        filename : string
            If set, a file will be written with all of the filtered halos
            and the quantities returned by the filter functions.
            Default=None.
        prefilters : array_like
            A single dataset can contain thousands or tens of thousands of
            halos. Significant time can be saved by not profiling halos
            that are certain to not pass any filter functions in place.
            Simple filters based on quantities provided in the initial
            halo list can be used to filter out unwanted halos using this
            parameter.
        
        Examples
        --------
        >>> hp.make_profiles(filename="FilteredQuantities.out",
                 prefilters=["halo['mass'] > 1e13"])
        
        """

        if len(self.all_halos) == 0:
            mylog.error("Halo list is empty, returning.")
            return None

        # Reset filtered halo list.
        self.filtered_halos = []

        # Check to see if the VirialFilter has been added to the filter list.
        # If a lower mass cutoff is being used, use it to make a pre-filter.
        if prefilters is None: prefilters = []
        virial_prefilter = None
        virial_prefilter_safety_factor = 0.5
        all_filter_functions = [hf['function'] for hf in self._halo_filters]
        virial_filter = VirialFilter in all_filter_functions
        if 'mass' in self.halo_list_format and VirialFilter in all_filter_functions:
            vFilter = self._halo_filters[all_filter_functions.index(VirialFilter)]
            if vFilter['kwargs'].has_key('virial_filters') and \
               vFilter['kwargs']['virial_filters'] is not None:
                all_vqFilters = [vqf[0] for vqf in vFilter['kwargs']['virial_filters']]
                if 'TotalMassMsun' in all_vqFilters:
                    mass_filter = vFilter['kwargs']['virial_filters'][all_vqFilters.index('TotalMassMsun')]
                    if '>' in mass_filter[1]:
                        virial_prefilter = "halo['mass'] %s %f * %s" % \
                            (mass_filter[1], virial_prefilter_safety_factor, mass_filter[2])
                        prefilters.append(virial_prefilter)
                    elif '<' in mass_filter[1]:
                        virial_prefilter = "halo['mass'] %s %f * %s" % \
                            (mass_filter[1], (1./virial_prefilter_safety_factor), mass_filter[2])
                        prefilters.append(virial_prefilter)

        # Add profile fields necessary for calculating virial quantities.
        if virial_filter: self._check_for_needed_profile_fields()

        # Create output directory.
        if self.output_dir is not None:
            self.__check_directory("%s/%s" % (self.output_dir, self.pf.directory))
            my_output_dir = "%s/%s/%s" % (self.output_dir, self.pf.directory, 
                                          self.profile_output_dir)
        else:
            my_output_dir = "%s/%s" % (self.pf.fullpath, self.profile_output_dir)
        self.__check_directory(my_output_dir)

        # Profile all halos.
        updated_halos = []
        for halo in parallel_objects(self.all_halos, -1):
            # Apply prefilters to avoid profiling unwanted halos.
            filter_result = True
            haloQuantities = {}
            if prefilters is not None:
                for prefilter in prefilters:
                    if not eval(prefilter):
                        filter_result = False
                        break

            if filter_result and len(self.profile_fields) > 0:

                profile_filename = "%s/Halo_%04d_profile.dat" % (my_output_dir, halo['id'])

                profiledHalo = self._get_halo_profile(halo, profile_filename, virial_filter=virial_filter)

                if profiledHalo is None:
                    continue

                # Apply filter and keep track of the quantities that are returned.
                for hFilter in self._halo_filters:
                    filter_result, filterQuantities = hFilter['function'](profiledHalo, *hFilter['args'], 
                                                                          **hFilter['kwargs'])

                    if not filter_result: break

                    if filterQuantities is not None:
                        haloQuantities.update(filterQuantities)

            if filter_result:
                for quantity in self.filter_quantities:
                    if halo.has_key(quantity): haloQuantities[quantity] = halo[quantity]

                self.filtered_halos.append(haloQuantities)

            # If we've gotten this far down, this halo is good and we want
            # to keep it. But we need to communicate the recentering changes
            # to all processors (the root one in particular) without having
            # one task clobber the other.
            updated_halos.append(halo)
        
        # And here is where we bring it all together.
        updated_halos = self.comm.par_combine_object(updated_halos,
                            datatype="list", op="cat")
        updated_halos.sort(key = lambda a:a['id'])
        self.all_halos = updated_halos

        self.filtered_halos = self.comm.par_combine_object(self.filtered_halos,
                            datatype="list", op="cat")
        self.filtered_halos.sort(key = lambda a:a['id'])

        if filename is not None:
            self._write_filtered_halo_list(filename, **kwargs)

    def _get_halo_profile(self, halo, filename, virial_filter=True,
            force_write=False):
        r"""Profile a single halo and write profile data to a file.
        If file already exists, read profile data from file.
        Return a dictionary of id, center, and virial quantities if virial_filter is True.
        """

        # Read profile from file if it already exists.
        # If not, profile will be None.
        profile = self._read_profile(filename)

        # Make profile if necessary.
        newProfile = profile is None
        if newProfile:

            r_min = 2 * self.pf.h.get_smallest_dx() * self.pf['mpc']
            if (halo['r_max'] / r_min < PROFILE_RADIUS_THRESHOLD):
                mylog.error("Skipping halo with r_max / r_min = %f." % (halo['r_max']/r_min))
                return None

            # get a sphere object to profile
            sphere = get_halo_sphere(halo, self.pf, recenter=self.recenter)
            if sphere is None: return None

            if self._need_bulk_velocity:
                # Set bulk velocity to zero out radial velocity profiles.
                if self.velocity_center[0] == 'bulk':
                    if self.velocity_center[1] == 'halo':
                        sphere.set_field_parameter('bulk_velocity', halo['velocity'])
                    elif self.velocity_center[1] == 'sphere':
                        sphere.set_field_parameter('bulk_velocity', 
                                                   sphere.quantities['BulkVelocity'](lazy_reader=False, 
                                                                                     preload=False))
                    else:
                        mylog.error("Invalid parameter: VelocityCenter.")
                elif self.velocity_center[0] == 'max':
                    max_grid, max_cell, max_value, max_location = \
                        self.pf.h.find_max_cell_location(self.velocity_center[1])
                    sphere.set_field_parameter('bulk_velocity', [max_grid['x-velocity'][max_cell],
                                                                 max_grid['y-velocity'][max_cell],
                                                                 max_grid['z-velocity'][max_cell]])

            try:
                profile = BinnedProfile1D(sphere, self.n_profile_bins, "RadiusMpc",
                                                r_min, halo['r_max'],
                                                log_space=True, lazy_reader=False,
                                                end_collect=True)
            except EmptyProfileData:
                mylog.error("Caught EmptyProfileData exception, returning None for this halo.")
                return None
            # Figure out which fields to add simultaneously
            field_groupings = defaultdict(lambda: defaultdict(list))
            for hp in self.profile_fields:
                field_groupings[hp['weight_field']][hp['accumulation']].append(hp['field'])
            for weight_field in field_groupings:
                for accum, fields in field_groupings[weight_field].items():
                    profile.add_fields(fields, weight=weight_field,
                                       accumulation=accum)

        if virial_filter:
            self._add_actual_overdensity(profile)

        if newProfile:
            mylog.info("Writing halo %d" % halo['id'])
            profile.write_out(filename, format='%0.6e')
        elif force_write:
            mylog.info("Re-writing halo %d" % halo['id'])
            self._write_profile(profile, filename, format='%0.6e')

        if newProfile:
            # Temporary solution to memory leak.
            for g in self.pf.h.grids:
                g.clear_data()
            sphere.clear_data()
            del sphere

        return profile

    @parallel_blocking_call
    def make_projections(self, axes=[0, 1, 2], halo_list='filtered',
            save_images=False, save_cube=True):
        r"""Make projections of all halos using specified fields.
        
        After adding fields using `add_projection`, this starts the actual
        calculations and saves the output to disk.
        
        Parameters
        ---------
        axes = array_like
            A list of the axes to project along, using the usual 0,1,2
            convention. Default=[0,1,2]
        halo_list : {'filtered', 'all'}
            Which set of halos to make profiles of, either ones passed by the
            halo filters (if enabled/added), or all halos.
            Default='filtered'.
        save_images : bool
            Whether or not to save images of the projections. Default=False.
        save_cube : bool
            Whether or not to save the HDF5 files of the halo projections.
            Default=True.
        
        Examples
        --------
        >>> hp.make_projections(axes=[0, 1, 2], save_cube=True,
            save_images=True, halo_list="filtered")
        
        """

        # Get list of halos for projecting.
        if halo_list == 'filtered':
            halo_projection_list = self.filtered_halos
        elif halo_list == 'all':
            halo_projection_list = self.all_halos
        elif isinstance(halo_list, types.StringType):
            halo_projection_list = self._read_halo_list(halo_list)
        elif isinstance(halo_list, types.ListType):
            halo_projection_list = halo_list
        else:
            mylog.error("Keyword, halo_list', must be 'filtered', 'all', a filename, or an actual list.")
            return

        if len(halo_projection_list) == 0:
            mylog.error("Halo list for projections is empty.")
            return

        # Set resolution for fixed resolution output.
        if self.project_at_level == 'max':
            proj_level = self.pf.h.max_level
        else:
            proj_level = int(self.project_at_level)
        proj_dx = self.pf.units[self.projection_width_units] / \
            self.pf.parameters['TopGridDimensions'][0] / \
            (self.pf.parameters['RefineBy']**proj_level)
        projectionResolution = int(self.projection_width / proj_dx)

        # Create output directory.
        if self.output_dir is not None:
            self.__check_directory("%s/%s" % (self.output_dir, self.pf.directory))
            my_output_dir = "%s/%s/%s" % (self.output_dir, self.pf.directory, 
                                          self.projection_output_dir)
        else:
            my_output_dir = "%s/%s" % (self.pf.fullpath, self.projection_output_dir)
        self.__check_directory(my_output_dir)

        center = [0.5 * (self.pf.parameters['DomainLeftEdge'][w] + 
                         self.pf.parameters['DomainRightEdge'][w])
                  for w in range(self.pf.parameters['TopGridRank'])]

        for halo in parallel_objects(halo_projection_list, -1):
            if halo is None:
                continue
            # Check if region will overlap domain edge.
            # Using non-periodic regions is faster than using periodic ones.
            leftEdge = [(halo['center'][w] - 
                         0.5 * self.projection_width/self.pf.units[self.projection_width_units])
                        for w in range(len(halo['center']))]
            rightEdge = [(halo['center'][w] + 
                          0.5 * self.projection_width/self.pf.units[self.projection_width_units])
                         for w in range(len(halo['center']))]

            mylog.info("Projecting halo %04d in region: [%f, %f, %f] to [%f, %f, %f]." %
                       (halo['id'], leftEdge[0], leftEdge[1], leftEdge[2], 
                        rightEdge[0], rightEdge[1], rightEdge[2]))

            need_per = False
            for w in range(len(halo['center'])):
                if ((leftEdge[w] < self.pf.parameters['DomainLeftEdge'][w]) or
                    (rightEdge[w] > self.pf.parameters['DomainRightEdge'][w])):
                    need_per = True
                    break

            if need_per:
                region = self.pf.h.periodic_region(halo['center'], leftEdge, rightEdge)
            else:
                region = self.pf.h.region(halo['center'], leftEdge, rightEdge)

            # Make projections.
            if not isinstance(axes, types.ListType): axes = list([axes])
            for w in axes:
                projections = []
                # YT projections do not follow the right-hand rule.
                coords = range(3)
                del coords[w]
                x_axis = coords[0]
                y_axis = coords[1]

                for hp in self.projection_fields:
                    projections.append(self.pf.h.proj(w, hp['field'], 
                                                      weight_field=hp['weight_field'], 
                                                      source=region, center=halo['center'],
                                                      serialize=False))
                
                # Set x and y limits, shift image if it overlaps domain boundary.
                if need_per:
                    pw = self.projection_width/self.pf.units[self.projection_width_units]
                    _shift_projections(self.pf, projections, halo['center'], center, w)
                    # Projection has now been shifted to center of box.
                    proj_left = [center[x_axis]-0.5*pw, center[y_axis]-0.5*pw]
                    proj_right = [center[x_axis]+0.5*pw, center[y_axis]+0.5*pw]
                else:
                    proj_left = [leftEdge[x_axis], leftEdge[y_axis]]
                    proj_right = [rightEdge[x_axis], rightEdge[y_axis]]

                # Save projection data to hdf5 file.
                if save_cube or save_images:
                    axis_labels = ['x', 'y', 'z']

                    if save_cube:
                        dataFilename = "%s/Halo_%04d_%s_data.h5" % \
                            (my_output_dir, halo['id'], axis_labels[w])
                        mylog.info("Saving projection data to %s." % dataFilename)
                        output = h5py.File(dataFilename, "a")

                    # Create fixed resolution buffer for each projection and write them out.
                    for e, hp in enumerate(self.projection_fields):
                        frb = FixedResolutionBuffer(projections[e], (proj_left[0], proj_right[0], 
                                                                     proj_left[1], proj_right[1]),
                                                    (projectionResolution, projectionResolution),
                                                    antialias=False)
                        dataset_name = "%s_%s" % (hp['field'], hp['weight_field'])
                        if save_cube:
                            if dataset_name in output: del output[dataset_name]
                            output.create_dataset(dataset_name, data=frb[hp['field']])
                        if save_images:
                            filename = "%s/Halo_%04d_%s_%s.png" % (my_output_dir, halo['id'], 
                                                                   dataset_name, axis_labels[w])
                            if (frb[hp['field']] != 0).any():
                                write_image(na.log10(frb[hp['field']]), filename, cmap_name=hp['cmap'])
                            else:
                                mylog.info('Projection of %s for halo %d is all zeros, skipping image.' %
                                            (hp['field'], halo['id']))
                    if save_cube: output.close()

            del region

    @parallel_blocking_call
    def analyze_halo_spheres(self, analysis_function, halo_list='filtered',
                             analysis_output_dir=None):
        r"""Perform custom analysis on all halos.
        
        This will loop through all halo on the HaloProfiler's list, 
        creating a sphere object for each halo and passing that sphere 
        to the provided analysis function.
        
        Parameters
        ---------
        analysis_function : function
            A function taking two arguments, the halo dictionary, and a 
            sphere object.
            Example function to calculate total mass of halo:
                def my_analysis(halo, sphere):
                    total_mass = sphere.quantities['TotalMass']()
                    print total_mass
        halo_list : {'filtered', 'all'}
            Which set of halos to make profiles of, either ones passed by the
            halo filters (if enabled/added), or all halos.
            Default='filtered'.
        analysis_output_dir : string, optional
            If specified, this directory will be created within the dataset to 
            contain any output from the analysis function.  Default: None.

        Examples
        --------
        >>> hp.analyze_halo_spheres(my_analysis, halo_list="filtered",
                                    analysis_output_dir='special_analysis')
        
        """

        # Get list of halos for projecting.
        if halo_list == 'filtered':
            halo_analysis_list = self.filtered_halos
        elif halo_list == 'all':
            halo_analysis_list = self.all_halos
        elif isinstance(halo_list, types.StringType):
            halo_analysis_list = self._read_halo_list(halo_list)
        elif isinstance(halo_list, types.ListType):
            halo_analysis_list = halo_list
        else:
            mylog.error("Keyword, halo_list', must be 'filtered', 'all', a filename, or an actual list.")
            return

        if len(halo_analysis_list) == 0:
            mylog.error("Halo list for analysis is empty.")
            return

        # Create output directory.
        if analysis_output_dir is not None:
            if self.output_dir is not None:
                self.__check_directory("%s/%s" % (self.output_dir, self.pf.directory))
                my_output_dir = "%s/%s/%s" % (self.output_dir, self.pf.directory, 
                                              analysis_output_dir)
            else:
                my_output_dir = "%s/%s" % (self.pf.fullpath, analysis_output_dir)
            self.__check_directory(my_output_dir)

        for halo in parallel_objects(halo_analysis_list, -1):
            if halo is None: continue

            # Get a sphere object to analze.
            sphere = get_halo_sphere(halo, self.pf, recenter=self.recenter)
            if sphere is None: continue

            # Call the given analysis function.
            analysis_function(halo, sphere)

    def _add_actual_overdensity(self, profile):
        "Calculate overdensity from TotalMassMsun and CellVolume fields."

        if 'ActualOverdensity' in profile.keys():
            return

        rho_crit_now = 1.8788e-29 * self.pf.hubble_constant**2 # g cm^-3
        Msun2g = 1.989e33
        rho_crit = rho_crit_now * ((1.0 + self.pf.current_redshift)**3.0)
        if not self.use_critical_density: rho_crit *= self.pf.omega_matter

        profile['ActualOverdensity'] = (Msun2g * profile['TotalMassMsun']) / \
            profile['CellVolume'] / rho_crit

    def _check_for_needed_profile_fields(self):
        "Make sure CellVolume and TotalMass fields are added so virial quantities can be calculated."
        all_profile_fields = [hp['field'] for hp in self.profile_fields]
        if not 'CellVolume' in all_profile_fields:
            mylog.info("Adding CellVolume field to so virial quantities can be calculated")
            self.add_profile('CellVolume', weight_field=None, accumulation=True)
        if not 'TotalMassMsun' in all_profile_fields:
            mylog.info("Adding TotalMassMsun field to so virial quantities can be calculated")
            self.add_profile('TotalMassMsun', weight_field=None, accumulation=True)

    def _load_halo_data(self, filename=None):
        "Read hop output file or run hop if it doesn't exist."

        # Don't run if hop data already loaded.
        if self.all_halos:
            return

        if filename is None:
            filename = self.halo_list_file

        if self.output_dir is not None:
            self.__check_directory("%s/%s" % (self.output_dir, self.pf.directory))
            hop_file = "%s/%s/%s" % (self.output_dir, self.pf.directory, filename)
        else:
            hop_file = "%s/%s" % (self.pf.fullpath, filename)

        if not(os.path.exists(hop_file)):
            mylog.info("Halo finder file not found, running halo finder to get halos.")
            self._run_hop(hop_file)

        self.all_halos = self._read_halo_list(hop_file)

    def _read_halo_list(self, listFile):
        """
        Read halo list from aue file.
        Allow for columnar data in varying formats.
        """

        def __isE(arg):
            parts = arg.lower().split('e')
            if len(parts) != 2: return False
            return not (True in [q.isalpha() for q in ''.join(parts)])

        def __get_num(arg):
            if __isE(arg):
                return float(arg)
            if arg != arg.swapcase():
                return arg
            return float(arg)

        mylog.info("Reading halo information from %s." % listFile)
        haloList = []
        listLines = file(listFile)

        fields = self.halo_list_format.keys()
        getID = not 'id' in fields
        has_rmax = 'r_max' in fields
        has_r200kpc = 'r200kpc' in fields

        for line in listLines:
            line = line.strip()
            if len(line) > 0 and not line.startswith('#') and not line[0].isalpha():
                halo = {}
                onLine = line.split()
                for field in fields:
                    if isinstance(self.halo_list_format[field], types.ListType):
                        halo[field] = [__get_num(onLine[q]) for q in self.halo_list_format[field]]
                    else:
                        halo[field] = __get_num(onLine[self.halo_list_format[field]])
                if getID: halo['id'] = len(haloList)
                if self.halo_radius is not None:
                    halo['r_max'] = self.halo_radius * self.pf.units['mpc']
                elif has_rmax:
                    halo['r_max'] *= self.pf.units['mpc']
                elif has_r200kpc:
                    # If P-Groupfinder used, r_200 [kpc] is calculated.
                    # set r_max as 50% past r_200.
                    halo['r_max'] = 1.5 * halo['r200kpc'] / 1000.
                else:
                    mylog.error("HaloProfiler has no way to get halo radius.")
                    return None
                haloList.append(halo)

        mylog.info("Loaded %d halos." % (len(haloList)))
        return haloList

    def _read_profile(self, profileFile):
        "Read radial profile from file.  Return None if it doesn't have all the fields requested."

        # Check to see if file exists.
        if not os.path.exists(profileFile):
            return None

        f = open(profileFile, 'r')
        lines = f.readlines()
        f.close()

        # Get fields from header.
        header = lines.pop(0)
        header = header.strip()
        fields = header.split()
        # First string is '#'.
        fields.pop(0)

        profile = {}
        profile_obj = FakeProfile(self.pf)
        for field in fields:
            profile[field] = []

        # Check if all fields needed are present.
        all_profile_fields = [hp['field'] for hp in self.profile_fields]
        for field in all_profile_fields:
            if not field in profile:
                return None

        # Fill profile fields, skip bad values.
        for line in lines:
            line = line.strip()
            onLine = line.split()
            lineOK = True
            for value in onLine:
                if value.isalpha():
                    lineOK = False
                    break
            if lineOK:
                for q, field in enumerate(fields):
                    profile[field].append(float(onLine[q]))

        for field in fields:
            profile[field] = na.array(profile[field])

        profile_obj._data = profile

        if len(profile[fields[0]]) > 1:
            return profile_obj
        else:
            return None

    @parallel_blocking_call
    def _run_hop(self, hop_file):
        "Run hop to get halos."

        hop_results = self.halo_finder_function(self.pf, *self.halo_finder_args, 
                                                **self.halo_finder_kwargs)
        hop_results.write_out(hop_file)

        del hop_results
        self.pf.h.clear_all_data()

    @parallel_root_only
    def _write_filtered_halo_list(self, filename, format="%s"):
        """
        Write out list of filtered halos along with any quantities 
        picked up during the filtering process.
        """

        if len(self.filtered_halos) == 0:
            mylog.error("No halos in filtered list.")
            return

        filename = "%s/%s" % (self.pf.fullpath, filename)
        mylog.info("Writing filtered halo list to %s." % filename)
        file = open(filename, "w")
        fields = [field for field in sorted(self.filtered_halos[0])]
        halo_fields = []
        for halo_field in self.filter_quantities:
            if halo_field in fields:
                fields.remove(halo_field)
                halo_fields.append(halo_field)
        # Make it so number of fields in header is same as number of data columns.
        header_fields = []
        for halo_field in halo_fields:
            if isinstance(self.filtered_halos[0][halo_field], types.ListType):
                header_fields.extend(["%s[%d]" % (halo_field, q) 
                                      for q in range(len(self.filtered_halos[0][halo_field]))])
            else:
                header_fields.append(halo_field)
        file.write("# ")
        file.write("\t".join(header_fields + fields + ["\n"]))

        for halo in self.filtered_halos:
            for halo_field in halo_fields:
                if isinstance(halo[halo_field], types.ListType):
                    field_data = na.array(halo[halo_field])
                    field_data.tofile(file, sep="\t", format=format)
                else:
                    if halo_field == 'id':
                        file.write("%04d" % halo[halo_field])
                    else:
                        file.write("%s" % halo[halo_field])
                file.write("\t")
            field_data = na.array([halo[field] for field in fields])
            field_data.tofile(file, sep="\t", format=format)
            file.write("\n")
        file.close()

    def _write_profile(self, profile, filename, format="%0.16e"):
        fid = open(filename, "w")
        fields = [field for field in sorted(profile.keys()) if field != "UsedBins"]
        fid.write("\t".join(["#"] + fields + ["\n"]))
        field_data = na.array([profile[field] for field in fields])
        for line in range(field_data.shape[1]):
            field_data[:, line].tofile(fid, sep="\t", format=format)
            fid.write("\n")
        fid.close()

    @parallel_root_only
    def __check_directory(self, my_output_dir):
        if (os.path.exists(my_output_dir)):
            if not(os.path.isdir(my_output_dir)):
                mylog.error("Output directory exists, but is not a directory: %s." % my_output_dir)
                raise IOError(my_output_dir)
        else:
            os.mkdir(my_output_dir)

def get_halo_sphere(halo, pf, recenter=None):
    r"""Returns a sphere object for a given halo.
        
    With a dictionary containing halo properties, such as center 
    and r_max, this creates a sphere object and optionally 
    recenters and recreates the sphere using a recentering function.
    This is to be used primarily to make spheres for a set of halos 
    loaded by the HaloProfiler.
    
    Parameters
    ----------
    halo : dict, required
        The dictionary containing halo properties used to make the sphere.
        Required entries:
            center : list with center coordinates.
            r_max : sphere radius in Mpc.
    pf : parameter file object, required
        The parameter file from which the sphere will be made.
    recenter : {None, string or function}
        The exact location of the sphere center can significantly affect 
        radial profiles.  The halo center loaded by the HaloProfiler will 
        typically be the dark matter center of mass calculated by a halo 
        finder.  However, this may not be the best location for centering 
        profiles of baryon quantities.  For example, one may want to center 
        on the maximum density.
        If recenter is given as a string, one of the existing recentering 
        functions will be used:
            Min_Dark_Matter_Density : location of minimum dark matter density
            Max_Dark_Matter_Density : location of maximum dark matter density
            CoM_Dark_Matter_Density : dark matter center of mass
            Min_Gas_Density : location of minimum gas density
            Max_Gas_Density : location of maximum gas density
            CoM_Gas_Density : gas center of mass
            Min_Total_Density : location of minimum total density
            Max_Total_Density : location of maximum total density
            CoM_Total_Density : total center of mass
            Min_Temperature : location of minimum temperature
            Max_Temperature : location of maximum temperature
        Alternately, a function can be supplied for custom recentering.
        The function should take only one argument, a sphere object.
            Example function:
                def my_center_of_mass(data):
                   my_x, my_y, my_z = data.quantities['CenterOfMass']()
                   return (my_x, my_y, my_z)

        Examples: this should primarily be used with the halo list of the HaloProfiler.
        This is an example with an abstract halo asssuming a pre-defined pf.
        >>> halo = {'center': [0.5, 0.5, 0.5], 'r_max': 1.0}
        >>> my_sphere = get_halo_sphere(halo, pf, recenter='Max_Gas_Density')
        >>> # Assuming the above example function has been defined.
        >>> my_sphere = get_halo_sphere(halo, pf, recenter=my_center_of_mass)
    """
        
    sphere = pf.h.sphere(halo['center'], halo['r_max']/pf.units['mpc'])
    if len(sphere._grids) == 0: return None
    new_sphere = False

    if recenter:
        old = halo['center']
        if recenter in centering_registry:
            new_x, new_y, new_z = \
                centering_registry[recenter](sphere)
        else:
            # user supplied function
            new_x, new_y, new_z = recenter(sphere)
        if new_x < pf.domain_left_edge[0] or \
                new_y < pf.domain_left_edge[1] or \
                new_z < pf.domain_left_edge[2]:
            mylog.info("Recentering rejected, skipping halo %d" % \
                halo['id'])
            return None
        halo['center'] = [new_x, new_y, new_z]
        d = pf['kpc'] * periodic_dist(old, halo['center'],
            pf.domain_right_edge - pf.domain_left_edge)
        mylog.info("Recentered halo %d %1.3e kpc away." % (halo['id'], d))
        # Expand the halo to account for recentering. 
        halo['r_max'] += d / 1000 # d is in kpc -> want mpc
        new_sphere = True

    if new_sphere:
        # Temporary solution to memory leak.
        for g in pf.h.grids:
            g.clear_data()
        sphere.clear_data()
        del sphere
        sphere = pf.h.sphere(halo['center'], halo['r_max']/pf.units['mpc'])
    return sphere

def _shift_projections(pf, projections, oldCenter, newCenter, axis):
    """
    Shift projection data around.
    This is necessary when projecting a preiodic region.
    """
    offset = [newCenter[q]-oldCenter[q] for q in range(len(oldCenter))]
    width = [pf.parameters['DomainRightEdge'][q]-pf.parameters['DomainLeftEdge'][q] \
                 for q in range(len(oldCenter))]

    del offset[axis]
    del width[axis]

    for plot in projections:
        # Get name of data field.
        other_fields = {'px':True, 'py':True, 'pdx':True, 'pdy':True, 'weight_field':True}
        for pfield in plot.field_data.keys():
            if not(other_fields.has_key(pfield)):
                field = pfield
                break

        # Shift x and y positions.
        plot['px'] += offset[0]
        plot['py'] += offset[1]

        # Wrap off-edge cells back around to other side (periodic boundary conditions).
        plot['px'][plot['px'] < 0] += width[0]
        plot['py'][plot['py'] < 0] += width[1]
        plot['px'][plot['px'] > width[0]] -= width[0]
        plot['py'][plot['py'] > width[1]] -= width[1]

        # After shifting, some cells have fractional coverage on both sides of the box.
        # Find those cells and make copies to be placed on the other side.

        # Cells hanging off the right edge.
        add_x_px = plot['px'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_px -= width[0]
        add_x_py = plot['py'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_pdx = plot['pdx'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_pdy = plot['pdy'][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_field = plot[field][plot['px'] + 0.5 * plot['pdx'] > width[0]]
        add_x_weight_field = plot['weight_field'][plot['px'] + 0.5 * plot['pdx'] > width[0]]

        # Cells hanging off the left edge.
        add2_x_px = plot['px'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_px += width[0]
        add2_x_py = plot['py'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_pdx = plot['pdx'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_pdy = plot['pdy'][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_field = plot[field][plot['px'] - 0.5 * plot['pdx'] < 0]
        add2_x_weight_field = plot['weight_field'][plot['px'] - 0.5 * plot['pdx'] < 0]

        # Cells hanging off the top edge.
        add_y_px = plot['px'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_py = plot['py'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_py -= width[1]
        add_y_pdx = plot['pdx'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_pdy = plot['pdy'][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_field = plot[field][plot['py'] + 0.5 * plot['pdy'] > width[1]]
        add_y_weight_field = plot['weight_field'][plot['py'] + 0.5 * plot['pdy'] > width[1]]

        # Cells hanging off the bottom edge.
        add2_y_px = plot['px'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_py = plot['py'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_py += width[1]
        add2_y_pdx = plot['pdx'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_pdy = plot['pdy'][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_field = plot[field][plot['py'] - 0.5 * plot['pdy'] < 0]
        add2_y_weight_field = plot['weight_field'][plot['py'] - 0.5 * plot['pdy'] < 0]

        # Add the hanging cells back to the projection data.
        plot.field_data['px'] = na.concatenate([plot['px'], add_x_px, add_y_px, 
                                                add2_x_px, add2_y_px])
        plot.field_data['py'] = na.concatenate([plot['py'], add_x_py, add_y_py, 
                                                add2_x_py, add2_y_py])
        plot.field_data['pdx'] = na.concatenate([plot['pdx'], add_x_pdx, add_y_pdx, 
                                                 add2_x_pdx, add2_y_pdx])
        plot.field_data['pdy'] = na.concatenate([plot['pdy'], add_x_pdy, add_y_pdy, 
                                                 add2_x_pdy, add2_y_pdy])
        plot.field_data[field] = na.concatenate([plot[field], add_x_field, add_y_field, 
                                                 add2_x_field, add2_y_field])
        plot.field_data['weight_field'] = na.concatenate([plot['weight_field'],
                                                          add_x_weight_field, add_y_weight_field, 
                                                          add2_x_weight_field, add2_y_weight_field])

        # Delete original copies of hanging cells.
        del add_x_px, add_y_px, add2_x_px, add2_y_px
        del add_x_py, add_y_py, add2_x_py, add2_y_py
        del add_x_pdx, add_y_pdx, add2_x_pdx, add2_y_pdx
        del add_x_pdy, add_y_pdy, add2_x_pdy, add2_y_pdy
        del add_x_field, add_y_field, add2_x_field, add2_y_field
        del add_x_weight_field, add_y_weight_field, add2_x_weight_field, add2_y_weight_field

class FakeProfile(ParallelAnalysisInterface):
    """
    This is used to mimic a profile object when reading profile data from disk.
    """
    def __init__(self, pf):
        ParallelAnalysisInterface.__init__(self)
        self.pf = pf
        self._data = {}

    def __getitem__(self, key):
        return self._data[key]

    def keys(self):
        return self._data.keys()

standard_fields = [
    ("Density", "CellMassMsun", False),
    ("Temperature", "CellMassMsun", False),
    ("VelocityMagnitude", "CellMassMsun", False),
    ("Ones", None, False),
    ("Entropy", "CellMassMsun", False),
    ("RadialVelocity", "CellMassMsun", False),
    ("SpecificAngularMomentumX", "CellMassMsun", False),
    ("SpecificAngularMomentumY", "CellMassMsun", False),
    ("SpecificAngularMomentumZ", "CellMassMsun", False),
    ("CoolingTime", "CellMassMsun", False),
    ("DynamicalTime", "CellMassMsun", False),
    ("CellMassMsun", None, True),
    ("TotalMassMsun", None, True),
    ("Dark_Matter_Density", "CellMassMsun", False),
    #("ParticleSpecificAngularMomentumX", "ParticleMassMsun"),
    #("ParticleSpecificAngularMomentumY", "ParticleMassMsun"),
    #("ParticleSpecificAngularMomentumZ", "ParticleMassMsun"),
    ("OverDensity", "CellMassMsun", False),
    #("ParticleMassMsun", None),
    ("StarParticleDensity", "StarParticleMassMsun", False), # How do we weight this?
    #("StarParticleMassMsun", None), 
    ("StarParticleDensity", "StarParticleMassMsun", False), # How do we weight this?
]

standard_fields += [("%s_Fraction" % (s), "CellMassMsun", False)
    for s in ["HI","HII","HeI","HeII","HeIII","H2I","H2II",
    "HM","Electron", "DI","DII","HDI","Metal"]
]


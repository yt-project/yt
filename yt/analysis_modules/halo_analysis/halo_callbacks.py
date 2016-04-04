"""
Halo callback object



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os

from yt.analysis_modules.cosmological_observation.light_ray.light_ray import \
    periodic_distance
from yt.data_objects.profiles import \
    create_profile
from yt.frontends.ytdata.utilities import \
    _hdf5_yt_array, \
    _yt_array_hdf5
from yt.units.yt_array import \
    YTArray
from yt.utilities.exceptions import \
    YTSphereTooSmall
from yt.funcs import \
    ensure_list
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.operator_registry import \
    OperatorRegistry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only
from yt.visualization.profile_plotter import \
    PhasePlot

callback_registry = OperatorRegistry()
    
def add_callback(name, function):
    callback_registry[name] =  HaloCallback(function)

class HaloCallback(object):
    r"""
    A HaloCallback is a function that minimally takes in a Halo object 
    and performs some analysis on it.  This function may attach attributes 
    to the Halo object, write out data, etc, but does not return anything.
    """
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, halo):
        self.function(halo, *self.args, **self.kwargs)
        return True

def halo_sphere(halo, radius_field="virial_radius", factor=1.0, 
                field_parameters=None):
    r"""
    Create a sphere data container to associate with a halo.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    radius_field : string
        Field to be retrieved from the quantities dictionary as 
        the basis of the halo radius.
        Default: "virial_radius".
    factor : float
        Factor to be multiplied by the base radius for defining 
        the radius of the sphere.
        Default: 1.0.
    field_parameters : dict
        Dictionary of field parameters to be set with the sphere 
        created.
        
    """

    dds = halo.halo_catalog.data_ds
    center = dds.arr([halo.quantities["particle_position_%s" % axis] \
                      for axis in "xyz"])
    radius = factor * halo.quantities[radius_field]
    if radius <= 0.0:
        halo.data_object = None
        return
    try:
        sphere = dds.sphere(center, radius)
    except YTSphereTooSmall:
        halo.data_object = None
        return
    if field_parameters is not None:
        for field, par in field_parameters.items():
            if isinstance(par, tuple) and par[0] == "quantity":
                value = halo.quantities[par[1]]
            else:
                value = par
            sphere.set_field_parameter(field, value)
    halo.data_object = sphere

add_callback("sphere", halo_sphere)

def sphere_field_max_recenter(halo, field):
    r"""
    Recenter the halo sphere on the location of the maximum of the given field.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    field : string
        Field to be used for recentering.

    """

    if halo.data_object is None: return
    s_ds = halo.data_object.ds
    old_sphere = halo.data_object
    max_vals = old_sphere.quantities.max_location(field)
    new_center = s_ds.arr(max_vals[1:])
    new_sphere = s_ds.sphere(new_center.in_units("code_length"),
                               old_sphere.radius.in_units("code_length"))
    mylog.info("Moving sphere center from %s to %s." % (old_sphere.center,
                                                        new_sphere.center))
    for par, value in old_sphere.field_parameters.items():
        if par not in new_sphere.field_parameters:
            new_sphere.set_field_parameter(par, value)
    halo.data_object = new_sphere

add_callback("sphere_field_max_recenter", sphere_field_max_recenter)

def sphere_bulk_velocity(halo):
    r"""
    Set the bulk velocity for the sphere.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.

    """

    halo.data_object.set_field_parameter("bulk_velocity",
        halo.data_object.quantities.bulk_velocity(use_particles=True))

add_callback("sphere_bulk_velocity", sphere_bulk_velocity)

def profile(halo, bin_fields, profile_fields, n_bins=32, extrema=None, logs=None, units=None,
            weight_field="cell_mass", accumulation=False, fractional=False,
            storage="profiles", output_dir="."):
    r"""
    Create 1, 2, or 3D profiles of a halo.

    Store profile data in a dictionary associated with the halo object.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    bin_fields : list of strings
        The binning fields for the profile.
    profile_fields : string or list of strings
        The fields to be profiled.
    n_bins : int or list of ints
        The number of bins in each dimension.  If None, 32 bins for
        each bin are used for each bin field.
        Default: 32.
    extrema : dict of min, max tuples
        Minimum and maximum values of the bin_fields for the profiles.
        The keys correspond to the field names. Defaults to the extrema
        of the bin_fields of the dataset. If a units dict is provided, extrema
        are understood to be in the units specified in the dictionary.
    logs : dict of boolean values
        Whether or not to log the bin_fields for the profiles.
        The keys correspond to the field names. Defaults to the take_log
        attribute of the field.
    units : dict of strings
        The units of the fields in the profiles, including the bin_fields.
    weight_field : string
        Weight field for profiling.
        Default : "cell_mass"
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of
        all the values from bin 0 to n.  If -True, the sum is reversed so
        that the value for bin n is the cumulative sum from bin N (total bins)
        to n.  If the profile is 2D or 3D, a list of values can be given to
        control the summation in each dimension independently.
        Default: False.
    fractional : If True the profile values are divided by the sum of all
        the profile data such that the profile represents a probability
        distribution function.
    storage : string
        Name of the dictionary to store profiles.
        Default: "profiles"
    output_dir : string
        Name of directory where profile data will be written.  The full path will be
        the output_dir of the halo catalog concatenated with this directory.
        Default : "."

    """

    mylog.info("Calculating 1D profile for halo %d." % 
               halo.quantities["particle_identifier"])

    dds = halo.halo_catalog.data_ds

    if dds is None:
        raise RuntimeError("Profile callback requires a data ds.")

    if not hasattr(halo, "data_object"):
        raise RuntimeError("Profile callback requires a data container.")

    if halo.data_object is None:
        mylog.info("Skipping halo %d since data_object is None." %
                   halo.quantities["particle_identifier"])
        return

    if output_dir is None:
        output_dir = storage
    output_dir = os.path.join(halo.halo_catalog.output_dir, output_dir)
    
    bin_fields = ensure_list(bin_fields)
    my_profile = create_profile(halo.data_object, bin_fields, profile_fields, n_bins=n_bins,
                                extrema=extrema, logs=logs, units=units, weight_field=weight_field,
                                accumulation=accumulation, fractional=fractional)
                  
    prof_store = dict([(field, my_profile[field]) \
                       for field in my_profile.field_data])
    prof_store[my_profile.x_field] = my_profile.x
    if len(bin_fields) > 1:
        prof_store[my_profile.y_field] = my_profile.y
    if len(bin_fields) > 2:
        prof_store[my_profile.z_field] = my_profile.z
    if hasattr(halo, storage):
        halo_store = getattr(halo, storage)
        if "used" in halo_store:
            halo_store["used"] &= my_profile.used
    else:
        halo_store = {"used": my_profile.used}
        setattr(halo, storage, halo_store)
    halo_store.update(prof_store)

    if hasattr(my_profile, "variance"):
        variance_store = dict([(field, my_profile.variance[field]) \
                           for field in my_profile.variance])
        variance_storage = "%s_variance" % storage
        if hasattr(halo, variance_storage):
            halo_variance_store = getattr(halo, variance_storage)
        else:
            halo_variance_store = {}
            setattr(halo, variance_storage, halo_variance_store)
        halo_variance_store.update(variance_store)

add_callback("profile", profile)

@parallel_root_only
def save_profiles(halo, storage="profiles", filename=None,
                  output_dir="."):
    r"""
    Save profile data to disk.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    storage : string
        Name of the dictionary attribute containing the profile data to be written.
        Default: "profiles"
    filename : string
        The name of the file to be written.  The final filename will be 
        "<filename>_<id>.h5".  If None, filename is set to the value given 
        by the storage keyword.
        Default: None
    output_dir : string
        Name of directory where profile data will be written.  The full path will be
        the output_dir of the halo catalog concatenated with this directory.
        Default : "."
    
    """

    if not hasattr(halo, storage):
        return
    
    if filename is None:
        filename = storage
    output_file = os.path.join(halo.halo_catalog.output_dir, output_dir,
                               "%s_%06d.h5" % (filename, 
                                               halo.quantities["particle_identifier"]))
    mylog.info("Saving halo %d profile data to %s." %
               (halo.quantities["particle_identifier"], output_file))

    fh = h5py.File(output_file, "w")
    my_profile = getattr(halo, storage)
    profile_group = fh.create_group("profiles")
    for field in my_profile:
        # Don't write code units because we might not know those later.
        if isinstance(my_profile[field], YTArray):
            my_profile[field].convert_to_cgs()
        _yt_array_hdf5(profile_group, str(field), my_profile[field])
    variance_storage = "%s_variance" % storage
    if hasattr(halo, variance_storage):
        my_profile = getattr(halo, variance_storage)
        variance_group = fh.create_group("variance")
        for field in my_profile:
            # Don't write code units because we might not know those later.
            if isinstance(my_profile[field], YTArray):
                my_profile[field].convert_to_cgs()
            _yt_array_hdf5(variance_group, str(field), my_profile[field])
    fh.close()

add_callback("save_profiles", save_profiles)

def load_profiles(halo, storage="profiles", fields=None,
                  filename=None, output_dir="."):
    r"""
    Load profile data from disk.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    storage : string
        Name of the dictionary attribute to store profile data.
        Default: "profiles"
    fields : string or list of strings
        The fields to be loaded.  If None, all fields present will be loaded.
        Default : None
    filename : string
        The name of the file to be loaded.  The final filename will be 
        "<filename>_<id>.h5".  If None, filename is set to the value given 
        by the storage keyword.
        Default: None
    output_dir : string
        Name of directory where profile data will be read.  The full path will be
        the output_dir of the halo catalog concatenated with this directory.
        Default : "."
    
    """

    if filename is None:
        filename = storage
    output_file = os.path.join(halo.halo_catalog.output_dir, output_dir,
                               "%s_%06d.h5" % (filename, 
                                               halo.quantities["particle_identifier"]))
    if not os.path.exists(output_file):
        raise RuntimeError("Profile file not found: %s." % output_file)
    mylog.info("Loading halo %d profile data from %s." %
               (halo.quantities["particle_identifier"], output_file))

    fh = h5py.File(output_file, "r")
    if fields is None:
        profile_fields = fh["profiles"].keys()
    else:
        profile_fields = fields
    my_profile = {}
    my_group = fh["profiles"]
    for field in profile_fields:
        if field not in my_group:
            raise RuntimeError("%s field not present in %s." % (field, output_file))
        my_profile[field] = _hdf5_yt_array(my_group, field,
                                           ds=halo.halo_catalog.halos_ds)
    setattr(halo, storage, my_profile)
    
    if "variance" in fh:
        my_variance = {}
        my_group = fh["variance"]
        if fields is None:
            profile_fields = my_group.keys()
        for field in profile_fields:
            if field not in my_group:
                raise RuntimeError("%s field not present in %s." % (field, output_file))
            my_variance[field] = _hdf5_yt_array(my_group, field,
                                                ds=halo.halo_catalog.halos_ds)
        setattr(halo, "%s_variance" % storage, my_variance)
        
    fh.close()

add_callback("load_profiles", load_profiles)

def virial_quantities(halo, fields, 
                      overdensity_field=("gas", "overdensity"),
                      critical_overdensity=200,
                      profile_storage="profiles"):
    r"""
    Calculate the value of the given fields at the virial radius defined at 
    the given critical density by interpolating from radial profiles.

    Parameters
    ----------    
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    fields : string or list of strings
        The fields whose virial values are to be calculated.
    overdensity_field : string or tuple of strings
        The field used as the overdensity from which interpolation is done to 
        calculate virial quantities.
        Default: ("gas", "overdensity")
    critical_overdensity : float
        The value of the overdensity at which to evaulate the virial quantities.  
        Overdensity is with respect to the critical density.
        Default: 200
    profile_storage : string
        Name of the halo attribute that holds the profiles to be used.
        Default: "profiles"
    
    """

    mylog.info("Calculating virial quantities for halo %d." %
               halo.quantities["particle_identifier"])

    fields = ensure_list(fields)
    fields = [halo.data_object._determine_fields(field)[0]
              for field in fields]
    
    dds = halo.halo_catalog.data_ds
    profile_data = getattr(halo, profile_storage)

    if overdensity_field not in profile_data:
      raise RuntimeError("virial_quantities callback requires profile of %s." %
                         str(overdensity_field))

    overdensity = profile_data[overdensity_field]
    dfilter = np.isfinite(overdensity) & profile_data["used"] & (overdensity > 0)

    v_fields = {}
    for field in fields:
        if isinstance(field, tuple):
            my_field = field[-1]
        else:
            my_field = field
        v_fields[field] = my_field
        v_field = "%s_%d" % (my_field, critical_overdensity)
        if v_field not in halo.halo_catalog.quantities:
            halo.halo_catalog.quantities.append(v_field)
    vquantities = dict([("%s_%d" % (v_fields[field], critical_overdensity),
                         dds.quan(0, profile_data[field].units)) \
                        for field in fields])

    if dfilter.sum() < 2:
        halo.quantities.update(vquantities)
        return

    # find interpolation index
    # require a negative slope, but not monotonicity
    vod = overdensity[dfilter].to_ndarray()
    if (vod > critical_overdensity).all():
        if vod[-1] < vod[-2]:
            index = -2
        else:
            halo.quantities.update(vquantities)
            return
    elif (vod < critical_overdensity).all():
        if vod[0] > vod[1]:
            index = 0
        else:
            halo.quantities.update(vquantities)
            return            
    else:
        # take first instance of downward intersection with critical value
        intersections = (vod[:-1] >= critical_overdensity) & \
            (vod[1:] < critical_overdensity)
        if not intersections.any():
            halo.quantities.update(vquantities)
            return            
        index = np.where(intersections)[0][0]

    for field in fields:
        v_prof = profile_data[field][dfilter].to_ndarray()
        slope = np.log(v_prof[index + 1] / v_prof[index]) / \
          np.log(vod[index + 1] / vod[index])
        value = dds.quan(np.exp(slope * np.log(critical_overdensity / 
                                               vod[index])) * v_prof[index],
                         profile_data[field].units).in_cgs()
        vquantities["%s_%d" % (v_fields[field], critical_overdensity)] = value

    halo.quantities.update(vquantities)

add_callback("virial_quantities", virial_quantities)

def phase_plot(halo, output_dir=".", phase_args=None, phase_kwargs=None):
    r"""
    Make a phase plot for the halo object.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    output_dir : string
        Name of directory where profile data will be written.  The full path will be
        the output_dir of the halo catalog concatenated with this directory.
        Default : "."
    phase_args : list
        List of arguments to be given to PhasePlot.
    phase_kwargs : dict
        Dictionary of keyword arguments to be given to PhasePlot.

    """

    if phase_args is None:
        phase_args = []
    if phase_kwargs is None:
        phase_kwargs = {}

    try:
        plot = PhasePlot(halo.data_object, *phase_args, **phase_kwargs)
        plot.save(os.path.join(halo.halo_catalog.output_dir, output_dir,
                               "halo_%06d" % halo.quantities["particle_identifier"]))
    except ValueError:
        return

add_callback("phase_plot", phase_plot)

def delete_attribute(halo, attribute):
    r"""
    Delete attribute from halo object.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    attribute : string
        The attribute to be deleted.

    """

    if hasattr(halo, attribute):
        delattr(halo, attribute)

add_callback("delete_attribute", delete_attribute)

def iterative_center_of_mass(halo, radius_field="virial_radius", inner_ratio=0.1, step_ratio=0.9,
                             units="pc"):
    r"""
    Adjust halo position by iteratively recalculating the center of mass while 
    decreasing the radius.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    radius_field : string
        The halo quantity to be used as the radius for the sphere.
        Default: "virial_radius"
    inner_ratio : float
        The ratio of the smallest sphere radius used for calculating the center of 
        mass to the initial radius.  The sphere radius is reduced and center of mass 
        recalculated until the sphere has reached this size.
        Default: 0.1
    step_ratio : float
        The multiplicative factor used to reduce the radius of the sphere after the 
        center of mass is calculated.
        Default: 0.9
    units : str
        The units for printing out the distance between the initial and final centers.
        Default : "pc"
        
    """
    if inner_ratio <= 0.0 or inner_ratio >= 1.0:
        raise RuntimeError("iterative_center_of_mass: inner_ratio must be between 0 and 1.")
    if step_ratio <= 0.0 or step_ratio >= 1.0:
        raise RuntimeError("iterative_center_of_mass: step_ratio must be between 0 and 1.")

    center_orig = halo.halo_catalog.data_ds.arr([halo.quantities["particle_position_%s" % axis]
                                                 for axis in "xyz"])
    sphere = halo.halo_catalog.data_ds.sphere(center_orig, halo.quantities[radius_field])

    while sphere.radius > inner_ratio * halo.quantities[radius_field]:
        new_center = sphere.quantities.center_of_mass(use_gas=True, use_particles=True)
        sphere = sphere.ds.sphere(new_center, step_ratio * sphere.radius)

    distance = periodic_distance(center_orig.in_units("code_length").to_ndarray(),
                                 new_center.in_units("code_length").to_ndarray())
    distance = halo.halo_catalog.data_ds.quan(distance, "code_length")
    mylog.info("Recentering halo %d %f %s away." %
               (halo.quantities["particle_identifier"],
                distance.in_units(units), units))

    for i, axis in enumerate("xyz"):
        halo.quantities["particle_position_%s" % axis] = sphere.center[i]
    del sphere
    
add_callback("iterative_center_of_mass", iterative_center_of_mass)

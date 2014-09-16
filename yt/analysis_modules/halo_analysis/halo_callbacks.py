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

import h5py
import numpy as np
import os

from yt.data_objects.profiles import \
     Profile1D
from yt.units.yt_array import \
     YTArray, YTQuantity
from yt.utilities.exceptions import \
     YTSphereTooSmall
from yt.funcs import \
     ensure_list, is_root
from yt.utilities.exceptions import YTUnitConversionError
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
        Defautl: 1.0.
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
    new_center = s_ds.arr(max_vals[2:])
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

def profile(halo, x_field, y_fields, x_bins=32, x_range=None, x_log=True,
            weight_field="cell_mass", accumulation=False, storage="profiles",
            output_dir="."):
    r"""
    Create 1d profiles.

    Store profile data in a dictionary associated with the halo object.

    Parameters
    ----------
    halo : Halo object
        The Halo object to be provided by the HaloCatalog.
    x_field : string
        The binning field for the profile.
    y_fields : string or list of strings
        The fields to be propython
        filed.
    x_bins : int
        The number of bins in the profile.
        Default: 32
    x_range : (float, float)
        The range of the x_field.  If None, the extrema are used.
        Default: None
    x_log : bool
        Flag for logarithmmic binning.
        Default: True
    weight_field : string
        Weight field for profiling.
        Default : "cell_mass"
    accumulation : bool
        If True, profile data is a cumulative sum.
        Default : False
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
    
    if x_range is None:
        x_range = list(halo.data_object.quantities.extrema(x_field, non_zero=True))

    my_profile = Profile1D(halo.data_object, x_field, x_bins, 
                           x_range[0], x_range[1], x_log, 
                           weight_field=weight_field)
    my_profile.add_fields(ensure_list(y_fields))

    # temporary fix since profiles do not have units at the moment
    for field in my_profile.field_data:
        my_profile.field_data[field] = dds.arr(my_profile[field],
                                               dds.field_info[field].units)

    # accumulate, if necessary
    if accumulation:
        used = my_profile.used
        for field in my_profile.field_data:
            if weight_field is None:
                my_profile.field_data[field][used] = \
                    np.cumsum(my_profile.field_data[field][used])
            else:
                my_weight = my_profile.weight
                my_profile.field_data[field][used] = \
                  np.cumsum(my_profile.field_data[field][used] * my_weight[used]) / \
                  np.cumsum(my_weight[used])
                  
    # create storage dictionary
    prof_store = dict([(field, my_profile[field]) \
                       for field in my_profile.field_data])
    prof_store[my_profile.x_field] = my_profile.x

    if hasattr(halo, storage):
        halo_store = getattr(halo, storage)
        if "used" in halo_store:
            halo_store["used"] &= my_profile.used
    else:
        halo_store = {"used": my_profile.used}
        setattr(halo, storage, halo_store)
    halo_store.update(prof_store)

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

    out_file = h5py.File(output_file, "w")
    my_profile = getattr(halo, storage)
    for field in my_profile:
        # Don't write code units because we might not know those later.
        if isinstance(my_profile[field], YTArray):
            my_profile[field].convert_to_cgs()
        dataset = out_file.create_dataset(str(field), data=my_profile[field])
        units = ""
        if isinstance(my_profile[field], YTArray):
            units = str(my_profile[field].units)
        dataset.attrs["units"] = units
    out_file.close()

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

    out_file = h5py.File(output_file, "r")
    my_profile = {}
    if fields is None:
        fields = out_file.keys()
    for field in fields:
        if field not in out_file:
            raise RuntimeError("%s field not present in %s." % (field, output_file))
        units = ""
        if "units" in out_file[field].attrs:
            units = out_file[field].attrs["units"]
        if units == "dimensionless": units = ""
        my_profile[field] = halo.halo_catalog.halos_ds.arr(out_file[field].value, units)
    out_file.close()
    setattr(halo, storage, my_profile)

add_callback("load_profiles", load_profiles)

def virial_quantities(halo, fields, critical_overdensity=200,
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
    critical_density : float
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
    
    dds = halo.halo_catalog.data_ds
    profile_data = getattr(halo, profile_storage)

    if ("gas", "overdensity") not in profile_data:
      raise RuntimeError('virial_quantities callback requires profile of ("gas", "overdensity").')

    overdensity = profile_data[("gas", "overdensity")]
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

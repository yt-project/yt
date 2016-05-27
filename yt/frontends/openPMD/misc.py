# Copyright (c) 2015, Axel Huebl, Remi Lehe
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

import h5py as h5
import numpy as np
import re
import string
import collections # for isinstance
import sys, getopt, os.path
from yt.utilities.logger import ytLogger as mylog

openPMD = "1.0.0"

ext_list = [["ED-PIC", np.uint32(1)]]

def open_file(file_name):
    try:
        f = h5.File(file_name, "r")
        return(f)
    except:
        raise

def get_attr(f, name):
    """
    Try to access the path `name` in the file `f`
    Return the corresponding attribute if it is present
    """
    if name in list(f.attrs.keys()):
        return(True, f.attrs[name])
    else:
        return(False, None)

def test_record(g, r):
    """
    Checks if a record is valid

    Parameters
    ----------
    g : h5py.Group
        The group the record resides in

    r : string
        The name of the record.

    Returns
    -------
    An array with 2 elements :
    - The first element is 1 if an error occured, and 0 otherwise
    - The second element is 0 if a warning arised, and 0 otherwise
    """
    regEx = re.compile("^\w+$") # Python3 only: re.ASCII
    if regEx.match(r):
        # test component names
        result_array = np.array([0,0])
        if not is_scalar_record(g[r]) :
            for component_name in g[r]:
                if not regEx.match(component_name):
                    mylog.warning("openPMD: Component %s of record %s is NOT" \
                    " named properly (a-Z0-9_)!" %(component_name, g[r].name) )
                    result_array += np.array([1,0])
    else:
        mylog.warning("openPMD: Record %s is NOT named properly (a-Z0-9_)!" \
              %(r.name) )
        result_array = np.array([1,0])

    return(result_array)

def test_key(f, v, request, name):
    """
    Checks whether a key is present. A key can either be
    a h5py.Group or a h5py.Dataset.
    Returns an error if the key if absent and requested
    Returns a warning if the key if absent and recommended

    Parameters
    ----------
    f : an h5py.File or h5py.Group object
        The object in which to find the key

    v : bool
        Verbose option

    request : string
        Either "required", "recommended" or "optional"

    name : string
        The name of the key within this File, Group or DataSet

    Returns
    -------
    An array with 2 elements :
    - The first element is 1 if an error occured, and 0 otherwise
    - The second element is 0 if a warning arised, and 0 otherwise
    """
    valid = (name in list(f.keys()))
    if valid:
        if v:
            mylog.info("Key %s (%s) exists in `%s`!" %(name, request, str(f.name) ) )
        result_array = np.array([0,0])
    else:
        if request == "required":
            mylog.warning("openPMD: Key %s (%s) does NOT exist in `%s`!" \
            %(name, request, str(f.name)) )
            result_array = np.array([1, 0])
        elif request == "recommended":
            mylog.info("openPMD: Key %s (%s) does NOT exist in `%s`!" \
            %(name, request, str(f.name)) )
            result_array = np.array([0, 1])
        elif request == "optional":
            if v:
                mylog.info("openPMD: Key %s (%s) does NOT exist in `%s`!"  \
            %(name, request, str(f.name)) )
            result_array = np.array([0, 0])
        else :
            raise ValueError("Unrecognized string for `request` : %s" %request)

    return(result_array)

def test_attr(f, v, request, name, is_type=None, type_format=None):
    """
    Checks whether an attribute is present.
    Returns an error if the attribute if absent and requested
    Returns a warning if the attribute if absent and recommanded

    Parameters
    ----------
    f : an h5py.File, h5py.Group or h5py.DataSet object
        The object in which to find the key

    v : bool
        Verbose option

    request : string
        Either "required", "recommended" or "optional

    name : string
        The name of the attribute within this File, Group or DataSet

    is_type : (numpy or python) data type
        The type of the attribute. Default is "arbitrary" for None.
        Can be a list of data types where at least one data type must match
        but this list can not be combined with type_format.

    type_format: (numpy or python) data type
        Used with is_type to specify numpy ndarray dtypes or a
        base np.string_ format regex. Can be a list of data types
        for ndarrays where at least one data type must match.

    Returns
    -------
    An array with 2 elements :
    - The first element is 1 if an error occured, and 0 otherwise
    - The second element is 0 if a warning arised, and 0 otherwise
    """
    valid, value = get_attr(f, name)
    if valid:
        if v:
            mylog.info("openPMD: Attribute %s (%s) exists in `%s`! Type = %s, Value = %s" \
            %(name, request, str(f.name), type(value), str(value)) )

        # test type
        if is_type is not None:
            if not type_format is None and not is_type is np.string_ and \
               not isinstance(type_format, collections.Iterable):
                type_format = [type_format]
                type_format_names = map(lambda x: x.__name__, type_format)
            if not is_type is None and not isinstance(is_type, collections.Iterable):
                is_type = [is_type]
            is_type_names = map(lambda x: x.__name__, is_type)
            # add for each type in is_type -> wrong, need to add this at the comparison level!
            if type(value) in is_type:
                # np.string_ format or general ndarray dtype text
                if type(value) is np.string_ and type_format is not None:
                    regEx = re.compile(type_format) # Python3 only: re.ASCII
                    if regEx.match(value.decode()) :
                        result_array = np.array([0,0])
                    else:
                        mylog.warning("openPMD: Attribute %s in `%s` does not satisfy " \
                              "format ('%s' should be in format '%s')!" \
                              %(name, str(f.name), value.decode(), type_format ) )
                        result_array = np.array([1,0])
                # ndarray dtypes
                elif type(value) is np.ndarray:
                    if value.dtype.type in type_format:
                        result_array = np.array([0,0])
                    elif type_format is None:
                        result_array = np.array([0,0])
                    else:
                        mylog.warning("openPMD: Attribute %s in `%s` is not of type " \
                              "ndarray of '%s' (is ndarray of '%s')!" \
                              %(name, str(f.name), type_format_names, \
                              value.dtype.type.__name__) )
                        result_array = np.array([1,0])
                else:
                    result_array = np.array([0,0])
            else:
                print(
                 "Error: Attribute %s in `%s` is not of type '%s' (is '%s')!" \
                 %(name, str(f.name), str(is_type_names), \
                  type(value).__name__) )
                result_array = np.array([1,0])
        else: # is_type is None (== arbitrary)
            result_array = np.array([0,0])
    else:
        if request == "required":
            mylog.warning("openPMD: Attribute %s (%s) does NOT exist in `%s`!" \
            %(name, request, str(f.name)) )
            result_array = np.array([1, 0])
        elif request == "recommended":
            mylog.info("openPMD: Attribute %s (%s) does NOT exist in `%s`!" \
            %(name, request, str(f.name)) )
            result_array = np.array([0, 1])
        elif request == "optional":
            if v:
                mylog.info("openPMD: Attribute %s (%s) does NOT exist in `%s`!"  \
            %(name, request, str(f.name)) )
            result_array = np.array([0, 0])
        else :
            raise ValueError("Unrecognized string for `request` : %s" %request)

    return(result_array)

def is_scalar_record(r):
    """
    Checks if a record is a scalar record or not.

    Parameters
    ----------
    r : an h5py.Group or h5py.Dataset object
        the record that shall be tested

    Returns
    -------
    bool : true if the record is a scalar record, false if the record
           is either a vector or an other type of tensor record
    """
    if type(r) is h5.Group :
        # now it could be either a vector/tensor record
        # or a scalar record with a constant component

        valid, value = get_attr(r, "value")
        # constant components require a "value" and a "shape" attribute
        if valid :
            return True
        else:
            return False
    else :
        return True

def test_component(c, v) :
    """
    Checks if a record component defines all required attributes.

    Parameters
    ----------
    c : an h5py.Group or h5py.Dataset object
        the record component that shall be tested

    v : bool
        Verbose option

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """
    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([0,0])

    if type(c) is h5.Group :
        # since this check tests components, this must be a constant
        # component: requires "value" and "shape" attributes
        result_array += test_attr(c, v, "required", "value") # type can be arbitrary
        result_array += test_attr(c, v, "required", "shape", np.ndarray, np.uint64)

    # default attributes for all components
    result_array += test_attr(c, v, "required", "unitSI", np.float64)

    return(result_array)


def check_root_attr(f, v, pic):
    """
    Scan the root of the file and make sure that all the attributes are present

    Parameters
    ----------
    f : an h5py.File object
        The HDF5 file in which to find the attribute

    v : bool
        Verbose option

    pic : bool
        Whether to check for the ED-PIC extension attributes

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """
    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([0,0])

    # STANDARD.md
    #   required
    result_array += test_attr(f, v, "required", "openPMD", np.string_, "^[0-9]+\.[0-9]+\.[0-9]+$")
    result_array += test_attr(f, v, "required", "openPMDextension", np.uint32)
    result_array += test_attr(f, v, "required", "basePath", np.string_, "^\/data\/\%T\/$")
    result_array += test_attr(f, v, "required", "meshesPath", np.string_)
    result_array += test_attr(f, v, "required", "particlesPath", np.string_)
    result_array += test_attr(f, v, "required", "iterationEncoding", np.string_, "^groupBased|fileBased$")
    result_array += test_attr(f, v, "required", "iterationFormat", np.string_)

    # groupBased iteration encoding needs to match basePath
    if result_array[0] == 0 :
        if f.attrs["iterationEncoding"].decode() == "groupBased" :
            if f.attrs["iterationFormat"].decode() != f.attrs["basePath"].decode() :
                mylog.warning("openPMD: for groupBased iterationEncoding the basePath "
                      "and iterationFormat must match!")
                result_array += np.array([1,0])

    #   recommended
    result_array += test_attr(f, v, "recommended", "author", np.string_)
    result_array += test_attr(f, v, "recommended", "software", np.string_)
    result_array += test_attr(f, v, "recommended",
                              "softwareVersion", np.string_)
    result_array += test_attr(f, v, "recommended", "date", np.string_,
      "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} [\+|-][0-9]{4}$")

    #   optional
    result_array += test_attr(f, v, "optional", "comment", np.string_)

    # Extension: ED-PIC
    if pic:
        valid, extensionIDs = get_attr(f, "openPMDextension")
        if valid:
            if (ext_list[0][1] & extensionIDs) != extensionIDs:
                mylog.warning("openPMD: ID=%s for extension `%s` not found in " \
                      "`openPMDextension` (is %s)!" \
                     %(ext_list[0][1], ext_list[0][0], extensionIDs) )
                result_array += np.array([1,0])

    return(result_array)


def check_iterations(f, v, pic) :
    """
    Scan all the iterations present in the file, checking both
    the meshes and the particles

    Parameters
    ----------
    f : an h5py.File object
        The HDF5 file in which to find the attribute

    v : bool
        Verbose option

    pic : bool
        Whether to check for the ED-PIC extension attributes

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """

    # Find all the iterations
    format_error = False
    try :
        list_iterations = list(f['/data/'].keys())
    except KeyError :
        format_error = True
    else :
        # Check that these iterations are indeed encoded as integers
        for iteration in list_iterations :
            for character in iteration : # go through the string
                if not (character in string.digits) :
                    format_error = True
    # Detect any error and interrupt execution if one is found
    if format_error == True :
        mylog.warning("openPMD: it seems that the path of the data within the HDF5 file "
              "is not of the form '/data/%T/', where %T corresponds to an "
              "actual integer.")
        return(np.array([1, 0]))
    else :
        mylog.info("openPMD: Found %d iteration(s)" % len(list_iterations) )

    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([ 0, 0])

    # Loop over the iterations and check the meshes and the particles
    for iteration in list_iterations :
        result_array += check_base_path(f, iteration, v, pic)
        # Go deeper only if there is no error at this point
        if result_array[0] == 0 :
            result_array += check_meshes(f, iteration, v, pic)
            result_array += check_particles(f, iteration, v, pic)

    return(result_array)

def check_base_path(f, iteration, v, pic):
    """
    Scan the base_path that corresponds to this iteration

    Parameters
    ----------
    f : an h5py.File object
        The HDF5 file in which to find the attribute

    iteration : string representing an integer
        The iteration at which to scan the meshes

    v : bool
        Verbose option

    pic : bool
        Whether to check for the ED-PIC extension attributes

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """
    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([ 0, 0])

    # Find the path to the data
    base_path = ("/data/%s/" % iteration).encode('ascii')
    bp = f[base_path]

    # Check for the attributes of the STANDARD.md
    result_array += test_attr(bp, v, "required", "time", [np.float32, np.float64])
    result_array += test_attr(bp, v, "required", "dt", [np.float32, np.float64])
    result_array += test_attr(bp, v, "required", "timeUnitSI", np.float64)

    return(result_array)

def check_meshes(f, iteration, v, pic):
    """
    Scan all the meshes corresponding to one iteration

    Parameters
    ----------
    f : an h5py.File object
        The HDF5 file in which to find the attribute

    iteration : string representing an integer
        The iteration at which to scan the meshes

    v : bool
        Verbose option

    pic : bool
        Whether to check for the ED-PIC extension attributes

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """
    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([ 0, 0])

    # Find the path to the data
    base_path = "/data/%s/" % iteration
    valid, meshes_path = get_attr(f, "meshesPath")
    if not valid :
        mylog.warning("openPMD: `meshesPath` is missing or malformed in '/'")
        return( np.array([1, 0]) )
    meshes_path = meshes_path.decode()

    if os.path.join( base_path, meshes_path) != ( base_path + meshes_path ):
        mylog.warning("openPMD: `basePath`+`meshesPath` seems to be malformed "
            "(is `basePath` absolute and ends on a `/` ?)")
        return( np.array([1, 0]) )
    else:
        full_meshes_path = (base_path + meshes_path).encode('ascii')
        # Find all the meshes
        try:
            list_meshes = list(f[full_meshes_path].keys())
        except KeyError:
            list_meshes = []
    print( "Iteration %s : found %d meshes"
        %( iteration, len(list_meshes) ) )

    # Check for the attributes of the STANDARD.md
    for field_name in list_meshes :
        field = f[full_meshes_path + field_name.encode('ascii')]

        result_array += test_record(f[full_meshes_path], field_name)

        # General attributes of the record
        result_array += test_attr(field, v, "required",
                                  "unitDimension", np.ndarray, np.float64)
        result_array += test_attr(field, v, "required",
                                  "timeOffset", [np.float32, np.float64])
        result_array += test_attr(field, v, "required",
                                  "gridSpacing", np.ndarray, [np.float32, np.float64])
        result_array += test_attr(field, v, "required",
                                  "gridGlobalOffset", np.ndarray, [np.float32, np.float64])
        result_array += test_attr(field, v, "required",
                                  "gridUnitSI", np.float64)
        result_array += test_attr(field, v, "required",
                                  "dataOrder", np.string_)
        result_array += test_attr(field, v, "required",
                                  "axisLabels", np.ndarray, np.string_)
        # Specific check for geometry
        geometry_test = test_attr(field, v, "required", "geometry", np.string_)
        result_array += geometry_test
        # geometryParameters is required when using thetaMode
        if geometry_test[0] == 0 and field.attrs["geometry"] == b"thetaMode" :
            result_array += test_attr(field, v, "required",
                                            "geometryParameters", np.string_)
        # otherwise it is optional
        else :
            result_array += test_attr(field, v, "optional",
                                            "geometryParameters", np.string_)

        # Attributes of the record's components
        if is_scalar_record(field) :   # If the record is a scalar field
            result_array += test_component(field, v)
            result_array += test_attr(field, v,
                                "required", "position", np.ndarray, [np.float32, np.float64])
        else:                          # If the record is a vector field
            # Loop over the components
            for component_name in list(field.keys()) :
                component = field[component_name]
                result_array += test_component(component, v)
                result_array += test_attr(component, v,
                                "required", "position", np.ndarray, [np.float32, np.float64])

    # Check for the attributes of the PIC extension,
    # if asked to do so by the user
    if pic:

        # Check the attributes associated with the field solver
        result_array += test_attr(f[full_meshes_path], v, "required",
                                  "fieldSolver", np.string_)
        valid, field_solver = get_attr(f[full_meshes_path], "fieldSolver")
        if (valid == True) and (field_solver in ["other", "GPSTD"]) :
            result_array += test_attr(f[full_meshes_path], v, "required",
                                      "fieldSolverParameters", np.string_)

        # Check for the attributes associated with the field boundaries
        result_array += test_attr(f[full_meshes_path], v, "required",
                                "fieldBoundary", np.ndarray, np.string_)
        valid, field_boundary = get_attr(f[full_meshes_path], "fieldBoundary")
        if (valid == True) and (np.any(field_boundary == b"other")) :
            result_array += test_attr(f[full_meshes_path], v, "required",
                        "fieldBoundaryParameters", np.ndarray, np.string_)

        # Check for the attributes associated with the field boundaries
        result_array += test_attr(f[full_meshes_path], v, "required",
                                "particleBoundary", np.ndarray, np.string_)
        valid, particle_boundary = get_attr(f[full_meshes_path], "particleBoundary")
        if (valid == True) and (np.any(particle_boundary == b"other")) :
            result_array += test_attr(f[full_meshes_path], v, "required",
                    "particleBoundaryParameters", np.ndarray, np.string_)

        # Check the attributes associated with the current smoothing
        result_array += test_attr(f[full_meshes_path], v, "required",
                                  "currentSmoothing", np.string_)
        valid, current_smoothing = get_attr(f[full_meshes_path], "currentSmoothing")
        if (valid == True) and (current_smoothing != b"none") :
            result_array += test_attr(f[full_meshes_path], v, "required",
                        "currentSmoothingParameters", np.string_)

        # Check the attributes associated with the charge conservation
        result_array += test_attr(f[full_meshes_path], v, "required",
                                  "chargeCorrection", np.string_)
        valid, charge_correction = get_attr(f[full_meshes_path], "chargeCorrection")
        if valid == True and charge_correction != b"none":
            result_array += test_attr(f[full_meshes_path], v, "required",
                        "chargeCorrectionParameters", np.string_)

        # Check for the attributes of each record
        for field_name in list_meshes :
            field = f[full_meshes_path + field_name.encode('ascii')]
            result_array + test_attr(field, v, "required",
                                     "fieldSmoothing", np.string_)
            valid, field_smoothing = get_attr(field, "fieldSmoothing")
            if (valid == True) and (field_smoothing != b"none") :
                result_array += test_attr(field,v, "required",
                                    "fieldSmoothingParameters", np.string_)
    return(result_array)


def check_particles(f, iteration, v, pic) :
    """
    Scan all the particle data corresponding to one iteration

    Parameters
    ----------
    f : an h5py.File object
        The HDF5 file in which to find the attribute

    iteration : string representing an integer
        The iteration at which to scan the particle data

    v : bool
        Verbose option

    pic : bool
        Whether to check for the ED-PIC extension attributes

    Returns
    -------
    An array with 2 elements :
    - The first element is the number of errors encountered
    - The second element is the number of warnings encountered
    """
    # Initialize the result array
    # First element : number of errors
    # Second element : number of warnings
    result_array = np.array([ 0, 0])

    # Find the path to the data
    base_path = ("/data/%s/" % iteration).encode('ascii')
    valid, particles_path = get_attr(f, "particlesPath")
    if os.path.join( base_path, particles_path) !=  \
        ( base_path + particles_path ) :
        mylog.warning("openPMD: `basePath`+`meshesPath` seems to be malformed "
            "(is `basePath` absolute and ends on a `/` ?)")
        return( np.array([1, 0]) )
    else:
        full_particle_path = base_path + particles_path
        # Find all the particle species
        try:
            list_species = list(f[full_particle_path].keys())
        except KeyError:
            list_species = []
    print( "Iteration %s : found %d particle species"
        %( iteration, len(list_species) ) )

    # Go through all the particle species
    for species_name in list_species :
        species = f[full_particle_path + species_name.encode('ascii')]

        # Check all records for this species
        for species_record_name in species :
            result_array += test_record(species, species_record_name)

        # Check the position record of the particles
        result_array += test_key(species, v, "required", "position")

        # Check the position offset record of the particles
        result_array += test_key(species, v, "required", "positionOffset")
        if result_array[0] == 0 :
            position_dimensions = len(species["position"].keys())
            positionOffset_dimensions = len(species["positionOffset"].keys())
            if position_dimensions != positionOffset_dimensions :
                mylog.warning("openPMD: `position` (ndim=%s) and `positionOffset` " \
                      "(ndim=%s) do not have the same dimensions in " \
                      "species `%s`!" \
                      %(str(position_dimensions), \
                        str(positionOffset_dimensions),
                        species.name) )
                result_array += np.array([ 1, 0])

        # Check the particlePatches record of the particles
        patch_test = test_key(species, v, "recommended", "particlePatches")
        result_array += patch_test
        if result_array[0] == 0 and patch_test[1] == 0 :
            result_array += test_key(species["particlePatches"], v, "required",
                                     "numParticles")
            result_array += test_key(species["particlePatches"], v, "required",
                                     "numParticlesOffset")
            result_array += test_key(species["particlePatches"], v, "required",
                                     "offset")
            result_array += test_key(species["particlePatches"], v, "required",
                                     "extent")
            if result_array[0] == 0 :
                offset = species["particlePatches"]["offset"]
                extent = species["particlePatches"]["extent"]
                # Attributes of the components
                for component_name in list(species["position"].keys()) :
                    result_array += test_key( offset, v, "required",
                                              component_name)
                    result_array += test_key( extent, v, "required",
                                              component_name)
                    if result_array[0] == 0 :
                        dset_offset = offset[component_name]
                        result_array += test_component(dset_offset, v)
                        dset_extent = extent[component_name]
                        result_array += test_component(dset_extent, v)

        # Check the records required by the PIC extension
        if pic :
            result_array += test_key(species, v, "required", "momentum")
            result_array += test_key(species, v, "required", "charge")
            result_array += test_key(species, v, "required", "mass")
            result_array += test_key(species, v, "required", "weighting")
            result_array += test_key(species, v, "optional", "boundElectrons")
            result_array += test_key(species, v, "optional", "protonNumber")
            result_array += test_key(species, v, "optional", "neutronNumber")

        # Check the attributes associated with the PIC extension
        if pic :
            result_array += test_attr(species, v, "required",
                                      "particleShape", [np.float32, np.float64])
            result_array += test_attr(species, v, "required",
                                      "currentDeposition", np.string_)
            result_array += test_attr(species, v, "required",
                                      "particlePush", np.string_)
            result_array += test_attr(species, v, "required",
                                      "particleInterpolation", np.string_)

            # Check for the attributes associated with the particle smoothing
            result_array += test_attr(species, v, "required",
                                      "particleSmoothing", np.string_)
            valid, particle_smoothing = get_attr(species, "particleSmoothing")
            if valid == True and particle_smoothing != b"none":
                result_array += test_attr(species, v, "required",
                                "particleSmoothingParameters", np.string_)

        # Check attributes of each record of the particle
        for record in list(species.keys()) :
            # all records (but particlePatches) require units
            if record != "particlePatches":
                result_array += test_attr(species[record], v,
                        "required", "unitDimension", np.ndarray, np.float64)
                time_type = f[base_path].attrs["time"].dtype.type
                result_array += test_attr(species[record], v, "required",
                                          "timeOffset", time_type)
                if pic :
                    result_array += test_attr(species[record], v, "required",
                                              "weightingPower", np.float64)
                    result_array += test_attr(species[record], v, "required",
                                              "macroWeighted", np.uint32)
                # Attributes of the components
                if is_scalar_record( species[record] ) : # Scalar record
                    dset = species[record]
                    result_array += test_component(dset, v)
                else : # Vector record
                    # Loop over the components
                    for component_name in list(species[record].keys()):
                        dset = species[ os.path.join(record, component_name) ]
                        result_array += test_component(dset, v)

    return(result_array)

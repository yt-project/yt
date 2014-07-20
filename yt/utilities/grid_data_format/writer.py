"""
Writing yt data to a GDF file.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import sys
import h5py
import numpy as np

from yt import __version__ as yt_version
from yt.utilities.exceptions import YTGDFAlreadyExists

def write_to_gdf(ds, gdf_path, data_author=None, data_comment=None,
                 dataset_units=None, particle_type_name="dark_matter",
                 clobber=False):
    """
    Write a dataset to the given path in the Grid Data Format.

    Parameters
    ----------
    ds : Dataset object
        The yt data to write out.
    gdf_path : string
        The path of the file to output.
    data_author : string, optional
        The name of the author who wrote the data. Default: None.
    data_comment : string, optional
        A descriptive comment. Default: None.
    dataset_units : dictionary, optional
        A dictionary of (value, unit) tuples to set the default units
        of the dataset. Keys can be:
            "length_unit"
            "time_unit"
            "mass_unit"
            "velocity_unit"
            "magnetic_unit"
        If not specified, these will carry over from the parent
        dataset.
    particle_type_name : string, optional
        The particle type of the particles in the dataset. Default: "dark_matter"
    clobber : boolean, optional
        Whether or not to clobber an already existing file. If False, attempting
        to overwrite an existing file will result in an exception.

    Examples
    --------
    >>> dataset_units = {"length_unit":(1.0,"Mpc"),
    ...                  "time_unit":(1.0,"Myr")}
    >>> write_to_gdf(ds, "clumps.h5", data_author="John ZuHone",
    ...              dataset_units=dataset_units,
    ...              data_comment="My Really Cool Dataset", clobber=True)
    """

    f = _create_new_gdf(ds, gdf_path, data_author, data_comment,
                        dataset_units=dataset_units,
                        particle_type_name=particle_type_name, clobber=clobber)

    # now add the fields one-by-one
    for field_name in ds.field_list:
        _write_field_to_gdf(ds, f, field_name, particle_type_name)

    # don't forget to close the file.
    f.close()


def save_field(ds, field_name, field_parameters=None):
    """
    Write a single field associated with the dataset ds to the
    backup file.

    Parameters
    ----------
    ds : Dataset object
        The yt dataset that the field is associated with.
    field_name : string
        The name of the field to save.
    field_parameters : dictionary
        A dictionary of field parameters to set.
    """

    if isinstance(field_name, tuple):
        field_name = field_name[1]
    field_obj = ds._get_field_info(field_name)
    if field_obj.particle_type:
        print("Saving particle fields currently not supported.")
        return

    backup_filename = ds.backup_filename
    if os.path.exists(backup_filename):
        # backup file already exists, open it
        f = h5py.File(backup_filename, "r+")
    else:
        # backup file does not exist, create it
        f = _create_new_gdf(ds, backup_filename, data_author=None,
                            data_comment=None,
                            particle_type_name="dark_matter")

    # now save the field
    _write_field_to_gdf(ds, f, field_name, particle_type_name="dark_matter",
                        field_parameters=field_parameters)

    # don't forget to close the file.
    f.close()


def _write_field_to_gdf(ds, fhandle, field_name, particle_type_name,
                        field_parameters=None):

    # add field info to field_types group
    g = fhandle["field_types"]
    # create the subgroup with the field's name
    if isinstance(field_name, tuple):
        field_name = field_name[1]
    fi = ds._get_field_info(field_name)
    try:
        sg = g.create_group(field_name)
    except ValueError:
        print "Error - File already contains field called " + field_name
        sys.exit(1)

    # grab the display name and units from the field info container.
    display_name = fi.display_name
    units = fi.units

    # check that they actually contain something...
    if display_name:
        sg.attrs["field_name"] = display_name
    else:
        sg.attrs["field_name"] = field_name
    if units:
        sg.attrs["field_units"] = units
    else:
        sg.attrs["field_units"] = "None"
    # @todo: is this always true?
    sg.attrs["staggering"] = 0

    # now add actual data, grid by grid
    g = fhandle["data"]
    for grid in ds.index.grids:

        # set field parameters, if specified
        if field_parameters is not None:
            for k, v in field_parameters.iteritems():
                grid.set_field_parameter(k, v)

        grid_group = g["grid_%010i" % (grid.id - grid._id_offset)]
        particles_group = grid_group["particles"]
        pt_group = particles_group[particle_type_name]
        # add the field data to the grid group
        # Check if this is a real field or particle data.
        grid.get_data(field_name)
        if fi.particle_type:  # particle data
            pt_group[field_name] = grid[field_name].in_units(units)
        else:  # a field
            grid_group[field_name] = grid[field_name].in_units(units)


def _create_new_gdf(ds, gdf_path, data_author=None, data_comment=None,
                    dataset_units=None, particle_type_name="dark_matter",
                    clobber=False):

    # Make sure we have the absolute path to the file first
    gdf_path = os.path.abspath(gdf_path)

    # Is the file already there? If so, are we allowing
    # clobbering?
    if os.path.exists(gdf_path) and not clobber:
        raise YTGDFAlreadyExists(gdf_path)

    ###
    # Create and open the file with h5py
    ###
    f = h5py.File(gdf_path, "w")

    ###
    # "gridded_data_format" group
    ###
    g = f.create_group("gridded_data_format")
    g.attrs["data_software"] = "yt"
    g.attrs["data_software_version"] = yt_version
    if data_author is not None:
        g.attrs["data_author"] = data_author
    if data_comment is not None:
        g.attrs["data_comment"] = data_comment

    ###
    # "simulation_parameters" group
    ###
    g = f.create_group("simulation_parameters")
    g.attrs["refine_by"] = ds.refine_by
    g.attrs["dimensionality"] = ds.dimensionality
    g.attrs["domain_dimensions"] = ds.domain_dimensions
    g.attrs["current_time"] = ds.current_time
    g.attrs["domain_left_edge"] = ds.domain_left_edge
    g.attrs["domain_right_edge"] = ds.domain_right_edge
    g.attrs["unique_identifier"] = ds.unique_identifier
    g.attrs["cosmological_simulation"] = ds.cosmological_simulation
    # @todo: Where is this in the yt API?
    g.attrs["num_ghost_zones"] = 0
    # @todo: Where is this in the yt API?
    g.attrs["field_ordering"] = 0
    # @todo: not yet supported by yt.
    g.attrs["boundary_conditions"] = np.array([0, 0, 0, 0, 0, 0], 'int32')

    if ds.cosmological_simulation:
        g.attrs["current_redshift"] = ds.current_redshift
        g.attrs["omega_matter"] = ds.omega_matter
        g.attrs["omega_lambda"] = ds.omega_lambda
        g.attrs["hubble_constant"] = ds.hubble_constant

    if dataset_units is None:
        dataset_units = {}

    g = f.create_group("dataset_units")
    for u in ["length","time","mass","velocity","magnetic"]:
        unit_name = u+"_unit"
        if unit_name in dataset_units:
            value, units = dataset_units[unit_name]
        else:
            attr = getattr(ds, unit_name)
            value = float(attr)
            units = str(attr.units)
        d = g.create_dataset(unit_name, data=value)
        d.attrs["unit"] = units

    ###
    # "field_types" group
    ###
    g = f.create_group("field_types")

    ###
    # "particle_types" group
    ###
    g = f.create_group("particle_types")

    # @todo: Particle type iterator
    sg = g.create_group(particle_type_name)
    sg["particle_type_name"] = particle_type_name

    ###
    # root datasets -- info about the grids
    ###
    f["grid_dimensions"] = ds.index.grid_dimensions
    f["grid_left_index"] = np.array(
        [grid.get_global_startindex() for grid in ds.index.grids]
    ).reshape(ds.index.grid_dimensions.shape[0], 3)
    f["grid_level"] = ds.index.grid_levels.flat
    # @todo: Fill with proper values
    f["grid_parent_id"] = -np.ones(ds.index.grid_dimensions.shape[0])
    f["grid_particle_count"] = ds.index.grid_particle_count

    ###
    # "data" group -- where we should spend the most time
    ###

    g = f.create_group("data")
    for grid in ds.index.grids:
        # add group for this grid
        grid_group = g.create_group("grid_%010i" % (grid.id - grid._id_offset))
        # add group for the particles on this grid
        particles_group = grid_group.create_group("particles")
        pt_group = particles_group.create_group(particle_type_name)

    return f

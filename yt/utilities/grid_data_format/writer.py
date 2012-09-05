"""
Writing yt data to a GDF file.

Authors: Casey W. Stark <caseywstark@gmail.com>
Affiliation: UC Berkeley

Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Casey W. Stark.  All Rights Reserved.

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

import os

import h5py
import numpy as np

from yt import __version__ as yt_version


def write_to_gdf(pf, gdf_path, data_author=None, data_comment=None,
                 particle_type_name="dark_matter"):
    """
    Write a parameter file to the given path in the Grid Data Format.

    Parameters
    ----------
    pf : StaticOutput object
        The yt data to write out.
    gdf_path : string
        The path of the file to output.

    """
    # Make sure we have the absolute path to the file first
    gdf_path = os.path.abspath(gdf_path)

    # Stupid check -- is the file already there?
    # @todo: make this a specific exception/error.
    if os.path.exists(gdf_path):
        raise IOError("A file already exists in the location: %s. Please provide a new one or remove that file." % gdf_path)

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
    g.attrs["refine_by"] = pf.refine_by
    g.attrs["dimensionality"] = pf.dimensionality
    g.attrs["domain_dimensions"] = pf.domain_dimensions
    g.attrs["current_time"] = pf.current_time
    g.attrs["domain_left_edge"] = pf.domain_left_edge
    g.attrs["domain_right_edge"] = pf.domain_right_edge
    g.attrs["unique_identifier"] = pf.unique_identifier
    g.attrs["cosmological_simulation"] = pf.cosmological_simulation
    # @todo: Where is this in the yt API?
    g.attrs["num_ghost_zones"] = 0
    # @todo: Where is this in the yt API?
    g.attrs["field_ordering"] = 0
    # @todo: not yet supported by yt.
    g.attrs["boundary_conditions"] = np.array([0, 0, 0, 0, 0, 0], 'int32')

    if pf.cosmological_simulation:
        g.attrs["current_redshift"] = pf.current_redshift
        g.attrs["omega_matter"] = pf.omega_matter
        g.attrs["omega_lambda"] = pf.omega_lambda
        g.attrs["hubble_constant"] = pf.hubble_constant

    ###
    # "field_types" group
    ###
    g = f.create_group("field_types")

    # Which field list should we iterate over?
    for field_name in pf.h.field_list:
        # create the subgroup with the field's name
        sg = g.create_group(field_name)

        # grab the display name and units from the field info container.
        display_name = pf.field_info[field_name].display_name
        units = pf.field_info[field_name].get_units()

        # check that they actually contain something...
        if display_name:
            sg.attrs["field_name"] = display_name
        else:
            sg.attrs["field_name"] = field_name
        if units:
            sg.attrs["field_units"] = units
        else:
            sg.attrs["field_units"] = "None"
        # @todo: the values must be in CGS already right?
        sg.attrs["field_to_cgs"] = 1.0
        # @todo: is this always true?
        sg.attrs["staggering"] = 0

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
    f["grid_dimensions"] = pf.h.grid_dimensions
    f["grid_left_index"] = np.array(
            [g.get_global_startindex() for g in pf.h.grids]
    ).reshape(pf.h.grid_dimensions.shape[0], 3)
    f["grid_level"] = pf.h.grid_levels
    # @todo: Fill with proper values
    f["grid_parent_id"] = -np.ones(pf.h.grid_dimensions.shape[0])
    f["grid_particle_count"] = pf.h.grid_particle_count

    ###
    # "data" group -- where we should spend the most time
    ###
    g = f.create_group("data")

    for grid in pf.h.grids:
        # add group for this grid

        grid_group = g.create_group("grid_%010i" % grid.id)
        # add group for the particles on this grid
        particles_group = grid_group.create_group("particles")
        pt_group = particles_group.create_group(particle_type_name)

        # add the field data to the grid group
        for field_name in pf.h.field_list:
            # Check if this is a real field or particle data.
            field_obj = pf.field_info[field_name]

            if field_obj.particle_type:  # particle data
                pt_group[field_name] = grid.get_data(field_name)
            else:  # a field
                grid_group[field_name] = grid.get_data(field_name)

    # don't forget to close the file.
    f.close()

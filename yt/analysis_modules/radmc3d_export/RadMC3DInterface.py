"""
Code to export from yt to RadMC3D



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.utilities.lib.write_array import \
    write_3D_array, write_3D_vector_array


class RadMC3DLayer:
    '''

    This class represents an AMR "layer" of the style described in
    the radmc3d manual. Unlike yt grids, layers may not have more
    than one parent, so level L grids will need to be split up
    if they straddle two or more level L - 1 grids.

    '''
    def __init__(self, level, parent, unique_id, LE, RE, dim):
        self.level = level
        self.parent = parent
        self.LeftEdge = LE
        self.RightEdge = RE
        self.ActiveDimensions = dim
        self.id = unique_id

    def get_overlap_with(self, grid):
        '''

        Returns the overlapping region between two Layers,
        or a layer and a grid. RE < LE means in any direction
        means no overlap.

        '''
        LE = np.maximum(self.LeftEdge,  grid.LeftEdge)
        RE = np.minimum(self.RightEdge, grid.RightEdge)
        return LE, RE

    def overlaps(self, grid):
        '''

        Returns whether or not this layer overlaps a given grid

        '''
        LE, RE = self.get_overlap_with(grid)
        if np.any(RE <= LE):
            return False
        else:
            return True


class RadMC3DWriter:
    '''

    This class provides a mechanism for writing out data files in a format
    readable by radmc3d. Currently, only the ASCII, "Layer" style file format
    is supported. For more information please see the radmc3d manual at:
    http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d

    Parameters
    ----------

    ds : `Dataset`
        This is the dataset object corresponding to the
        simulation output to be written out.

    max_level : int
        An int corresponding to the maximum number of levels of refinement
        to include in the output. Often, this does not need to be very large
        as information on very high levels is frequently unobservable.
        Default = 2.

    Examples
    --------

    This will create a field called "DustDensity" and write it out to the
    file "dust_density.inp" in a form readable by RadMC3D.

    >>> import yt
    >>> from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter

    >>> dust_to_gas = 0.01
    >>> def _DustDensity(field, data):
    ...     return dust_to_gas*data["Density"]
    >>> yt.add_field("DustDensity", function=_DustDensity)

    >>> ds = yt.load("galaxy0030/galaxy0030")

    >>> writer = RadMC3DWriter(ds)
    >>> writer.write_amr_grid()
    >>> writer.write_dust_file("DustDensity", "dust_density.inp")

    ---

    This example will create a field called "NumberDensityCO" and write it out
    to the file "numberdens_co.inp". It will also write out information about
    the gas velocity to "gas_velocity.inp" so that this broadening may be
    included in the radiative transfer calculation by radmc3d:

    >>> import yt
    >>> from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter

    >>> x_co = 1.0e-4
    >>> mu_h = yt.Quantity(2.34e-24, 'g')
    >>> def _NumberDensityCO(field, data):
    ...     return (x_co/mu_h)*data["Density"]
    >>> yt.add_field("NumberDensityCO", function=_NumberDensityCO)

    >>> ds = yt.load("galaxy0030/galaxy0030")
    >>> writer = RadMC3DWriter(ds)

    >>> writer.write_amr_grid()
    >>> writer.write_line_file("NumberDensityCO", "numberdens_co.inp")
    >>> velocity_fields = ["velocity_x", "velocity_y", "velocity_z"]
    >>> writer.write_line_file(velocity_fields, "gas_velocity.inp")

    '''

    def __init__(self, ds, max_level=2):
        self.max_level = max_level
        self.cell_count = 0
        self.layers = []
        self.domain_dimensions = ds.domain_dimensions
        self.domain_left_edge = ds.domain_left_edge
        self.domain_right_edge = ds.domain_right_edge
        self.grid_filename = "amr_grid.inp"
        self.ds = ds

        base_layer = RadMC3DLayer(0, None, 0,
                                  self.domain_left_edge,
                                  self.domain_right_edge,
                                  self.domain_dimensions)

        self.layers.append(base_layer)
        self.cell_count += np.product(ds.domain_dimensions)

        sorted_grids = sorted(ds.index.grids, key=lambda x: x.Level)
        for grid in sorted_grids:
            if grid.Level <= self.max_level:
                self._add_grid_to_layers(grid)

    def _get_parents(self, grid):
        parents = []
        for potential_parent in self.layers:
            if potential_parent.level == grid.Level - 1:
                if potential_parent.overlaps(grid):
                    parents.append(potential_parent)
        return parents

    def _add_grid_to_layers(self, grid):
        parents = self._get_parents(grid)
        for parent in parents:
            LE, RE = parent.get_overlap_with(grid)
            N = (RE - LE) / grid.dds
            N = np.array([int(n + 0.5) for n in N])
            new_layer = RadMC3DLayer(grid.Level, parent.id,
                                     len(self.layers),
                                     LE, RE, N)
            self.layers.append(new_layer)
            self.cell_count += np.product(N)

    def write_amr_grid(self):
        '''
        This routine writes the "amr_grid.inp" file that describes the mesh
        radmc3d will use.

        '''
        dims = self.domain_dimensions
        LE = self.domain_left_edge
        RE = self.domain_right_edge

        # RadMC-3D wants the cell wall positions in cgs. Convert here:
        LE_cgs = LE.in_units('cm').d  # don't write the units, though
        RE_cgs = RE.in_units('cm').d

        # calculate cell wall positions
        xs = [str(x) for x in np.linspace(LE_cgs[0], RE_cgs[0], dims[0]+1)]
        ys = [str(y) for y in np.linspace(LE_cgs[1], RE_cgs[1], dims[1]+1)]
        zs = [str(z) for z in np.linspace(LE_cgs[2], RE_cgs[2], dims[2]+1)]

        # writer file header
        grid_file = open(self.grid_filename, 'w')
        grid_file.write('1 \n')  # iformat is always 1
        if self.max_level == 0:
            grid_file.write('0 \n')
        else:
            grid_file.write('10 \n')  # only layer-style files are supported
        grid_file.write('1 \n')  # only cartesian coordinates are supported
        grid_file.write('0 \n')
        grid_file.write('{}    {}    {} \n'.format(1, 1, 1))  # assume 3D
        grid_file.write('{}    {}    {} \n'.format(dims[0], dims[1], dims[2]))
        if self.max_level != 0:
            s = str(self.max_level) + '    ' + str(len(self.layers)-1) + '\n'
            grid_file.write(s)

        # write base grid cell wall positions
        for x in xs:
            grid_file.write(x + '    ')
        grid_file.write('\n')

        for y in ys:
            grid_file.write(y + '    ')
        grid_file.write('\n')

        for z in zs:
            grid_file.write(z + '    ')
        grid_file.write('\n')

        # write information about fine layers, skipping the base layer:
        for layer in self.layers[1:]:
            p = layer.parent
            dds = (layer.RightEdge - layer.LeftEdge) / (layer.ActiveDimensions)
            if p == 0:
                ind = (layer.LeftEdge - LE) / (2.0*dds) + 1
            else:
                parent_LE = np.zeros(3)
                for potential_parent in self.layers:
                    if potential_parent.id == p:
                        parent_LE = potential_parent.LeftEdge
                ind = (layer.LeftEdge - parent_LE) / (2.0*dds) + 1
            ix = int(ind[0]+0.5)
            iy = int(ind[1]+0.5)
            iz = int(ind[2]+0.5)
            nx, ny, nz = layer.ActiveDimensions / 2
            s = '{}    {}    {}    {}    {}    {}    {} \n'
            s = s.format(p, ix, iy, iz, nx, ny, nz)
            grid_file.write(s)

        grid_file.close()

    def _write_layer_data_to_file(self, fhandle, field, level, LE, dim):
        cg = self.ds.covering_grid(level, LE, dim, num_ghost_zones=1)
        if isinstance(field, list):
            data_x = cg[field[0]]
            data_y = cg[field[1]]
            data_z = cg[field[2]]
            write_3D_vector_array(data_x, data_y, data_z, fhandle)
        else:
            data = cg[field]
            write_3D_array(data, fhandle)

    def write_dust_file(self, field, filename):
        '''
        This method writes out fields in the format radmc3d needs to compute
        thermal dust emission. In particular, if you have a field called
        "DustDensity", you can write out a dust_density.inp file.

        Parameters
        ----------

        field : string
            The name of the field to be written out
        filename : string
            The name of the file to write the data to. The filenames radmc3d
            expects for its various modes of operations are described in the
            radmc3d manual.

        '''
        fhandle = open(filename, 'w')

        # write header
        fhandle.write('1 \n')
        fhandle.write(str(self.cell_count) + ' \n')
        fhandle.write('1 \n')

        # now write fine layers:
        for layer in self.layers:
            lev = layer.level
            if lev == 0:
                LE = self.domain_left_edge
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)

        fhandle.close()

    def write_line_file(self, field, filename):
        '''
        This method writes out fields in the format radmc3d needs to compute
        line emission.

        Parameters
        ----------

        field : string or list of 3 strings
            If a string, the name of the field to be written out. If a list,
            three fields that will be written to the file as a vector quantity.
        filename : string
            The name of the file to write the data to. The filenames radmc3d
            expects for its various modes of operation are described in the
            radmc3d manual.

        '''
        fhandle = open(filename, 'w')

        # write header
        fhandle.write('1 \n')
        fhandle.write(str(self.cell_count) + ' \n')

        # now write fine layers:
        for layer in self.layers:
            lev = layer.level
            if lev == 0:
                LE = self.domain_left_edge
                N = self.domain_dimensions
            else:
                LE = layer.LeftEdge
                N = layer.ActiveDimensions

            self._write_layer_data_to_file(fhandle, field, lev, LE, N)

        fhandle.close()

    def write_source_files(self, sources, wavelengths):
        '''

        This function creates the stars.inp and wavelength_micron.inp
        files that RadMC3D uses for its dust continuum calculations.

        Parameters
        ----------

        sources: a list of RadMC3DSource objects
            A list that contains all the sources you would like yt
            to create
        wavelengths: np.array of float values
            An array listing the wavelength points you would like to
            use the radiative transfer calculation

        '''

        nstars = len(sources)
        nlam = len(wavelengths)

        filename = 'stars.inp'
        fhandle = open(filename, 'w')

        # write header
        fhandle.write('2 \n')  # a format flag that should always be 2
        fhandle.write('%d    %d \n' % (nstars, nlam))

        # write source information
        for source in sources:
            fhandle.write(str(source.radius) + ' ')
            fhandle.write(str(source.mass) + ' ')
            fhandle.write('%f %f %f' %(source.position[0], \
                                       source.position[1], \
                                       source.position[2]))
            fhandle.write('\n')

        # write wavelength informaton
        for wavelength in wavelengths:
            fhandle.write('%f \n' % wavelength)

        # finally write blackbody temperature for each source
        for source in sources:
            # the negative sign is a flag used internally
            # by RadMC3D to indicate that this is a blackbody
            # source
            fhandle.write('%f \n' % -source.temperature)

        # done with stars.inp
        fhandle.close()

        # now do the wavelength_micron.inp file
        filename = 'wavelength_micron.inp'
        fhandle = open(filename, 'w')

        fhandle.write('%d \n' % nlam)
        for wavelength in wavelengths:
            fhandle.write('%f \n' % wavelength)

        # done with both
        fhandle.close()


class RadMC3DSource:
    '''

    A class that contains the data associated with a single RadMC3D photon source.
    This is designed to help export data about the stars in a dataset into a format 
    that can be read in by RadMC3D. Although RadMC3D can handle non-blackbody 
    sources, here we assume that the source is a blackbody with a given temperature.

    Parameters
    ----------

    radius: float
        The size of the source in cm
    mass: float 
        The mass of the source in g
    position: list of floats
        The x, y, and z coordinates of the source, in cm
    temperature: float
        The blackbody temperature of the source, in K

    '''

    def __init__(self, radius, mass, position, temperature):
        self.radius = radius
        self.mass = mass
        self.position = position
        self.temperature = temperature

        # some basic sanity checks
        assert(self.radius > 0.0)
        assert(self.mass > 0.0)
        assert(self.temperature > 0)
        assert(len(self.position) == 3)  # 3D only, please

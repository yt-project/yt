"""
Generating PPV FITS cubes
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.utilities.on_demand_imports import _astropy
from yt.utilities.orientation import Orientation
from yt.utilities.fits_image import FITSImageData, sanitize_fits_unit
from yt.visualization.volume_rendering.off_axis_projection import off_axis_projection
from yt.funcs import get_pbar
from yt.utilities.physical_constants import clight, mh
import yt.units.dimensions as ytdims
from yt.units.yt_array import YTQuantity
from yt.funcs import iterable
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only, parallel_objects
import re
from . import ppv_utils
from yt.funcs import is_root
from yt.extern.six import string_types

def create_vlos(normal, no_shifting):
    if no_shifting:
        def _v_los(field, data):
            return data.ds.arr(data["zeros"], "cm/s")
    elif isinstance(normal, string_types):
        def _v_los(field, data):
            return -data["velocity_%s" % normal]
    else:
        orient = Orientation(normal)
        los_vec = orient.unit_vectors[2]
        def _v_los(field, data):
            vz = data["velocity_x"]*los_vec[0] + \
                data["velocity_y"]*los_vec[1] + \
                data["velocity_z"]*los_vec[2]
            return -vz
    return _v_los

fits_info = {"velocity":("m/s","VOPT","v"),
             "frequency":("Hz","FREQ","f"),
             "energy":("eV","ENER","E"),
             "wavelength":("angstrom","WAVE","lambda")}

class PPVCube(object):
    def __init__(self, ds, normal, field, velocity_bounds, center="c",
                 width=(1.0,"unitary"), dims=100, thermal_broad=False,
                 atomic_weight=56., depth=(1.0,"unitary"), depth_res=256,
                 method="integrate", weight_field=None, no_shifting=False,
                 north_vector=None, no_ghost=True, data_source=None):
        r""" Initialize a PPVCube object.

        Parameters
        ----------
        ds : dataset
            The dataset.
        normal : array_like or string
            The normal vector along with to make the projections. If an array, it
            will be normalized. If a string, it will be assumed to be along one of the
            principal axes of the domain ("x", "y", or "z").
        field : string
            The field to project.
        velocity_bounds : tuple
            A 4-tuple of (vmin, vmax, nbins, units) for the velocity bounds to
            integrate over.
        center : A sequence of floats, a string, or a tuple.
            The coordinate of the center of the image. If set to 'c', 'center' or
            left blank, the plot is centered on the middle of the domain. If set to
            'max' or 'm', the center will be located at the maximum of the
            ('gas', 'density') field. Centering on the max or min of a specific
            field is supported by providing a tuple such as ("min","temperature") or
            ("max","dark_matter_density"). Units can be specified by passing in *center*
            as a tuple containing a coordinate and string unit name or by passing
            in a YTArray. If a list or unitless array is supplied, code units are
            assumed.
        width : float, tuple, or YTQuantity.
            The width of the projection. A float will assume the width is in code units.
            A (value, unit) tuple or YTQuantity allows for the units of the width to be
            specified. Implies width = height, e.g. the aspect ratio of the PPVCube's
            spatial dimensions is 1.
        dims : integer, optional
            The spatial resolution of the cube. Implies nx = ny, e.g. the
            aspect ratio of the PPVCube's spatial dimensions is 1.
        thermal_broad : boolean, optional
            Whether or not to broaden the line using the gas temperature. Default: False.
        atomic_weight : float, optional
            Set this value to the atomic weight of the particle that is emitting the line
            if *thermal_broad* is True. Defaults to 56 (Fe).
        depth : A tuple or a float, optional
            A tuple containing the depth to project through and the string
            key of the unit: (width, 'unit').  If set to a float, code units
            are assumed. Only for off-axis cubes.
        depth_res : integer, optional
            Deprecated, this is still in the function signature for API
            compatibility
        method : string, optional
            Set the projection method to be used.
            "integrate" : line of sight integration over the line element.
            "sum" : straight summation over the line of sight.
        weight_field : string, optional
            The name of the weighting field.  Set to None for no weight.
        no_shifting : boolean, optional
            If set, no shifting due to velocity will occur but only thermal broadening.
            Should not be set when *thermal_broad* is False, otherwise nothing happens!
        north_vector : a sequence of floats
            A vector defining the 'up' direction. This option sets the orientation of
            the plane of projection. If not set, an arbitrary grid-aligned north_vector
            is chosen. Ignored in the case of on-axis cubes.
        no_ghost: bool, optional
            Optimization option for off-axis cases. If True, homogenized bricks will
            extrapolate out from grid instead of interpolating from
            ghost zones that have to first be calculated.  This can
            lead to large speed improvements, but at a loss of
            accuracy/smoothness in resulting image.  The effects are
            less notable when the transfer function is smooth and
            broad. Default: True
        data_source : yt.data_objects.data_containers.YTSelectionContainer, optional
            If specified, this will be the data source used for selecting regions to project.

        Examples
        --------
        >>> i = 60*np.pi/180.
        >>> L = [0.0,np.sin(i),np.cos(i)]
        >>> cube = PPVCube(ds, L, "density", (-5.,4.,100,"km/s"), width=(10.,"kpc"))
        """

        self.ds = ds
        self.field = field
        self.width = width
        self.particle_mass = atomic_weight*mh
        self.thermal_broad = thermal_broad
        self.no_shifting = no_shifting

        if not isinstance(normal, string_types):
            width = ds.coordinates.sanitize_width(normal, width, depth)
            width = tuple(el.in_units('code_length').v for el in width)

        if no_shifting and not thermal_broad:
            raise RuntimeError("no_shifting cannot be True when thermal_broad is False!")

        self.center = ds.coordinates.sanitize_center(center, normal)[0]

        self.nx = dims
        self.ny = dims
        self.nv = velocity_bounds[2]

        if method not in ["integrate","sum"]:
            raise RuntimeError("Only the 'integrate' and 'sum' projection +"
                               "methods are supported in PPVCube.")

        dd = ds.all_data()
        fd = dd._determine_fields(field)[0]
        self.field_units = ds._get_field_info(fd).units

        self.vbins = ds.arr(np.linspace(velocity_bounds[0],
                                        velocity_bounds[1],
                                        velocity_bounds[2]+1), velocity_bounds[3])

        self._vbins = self.vbins.copy()
        self.vmid = 0.5*(self.vbins[1:]+self.vbins[:-1])
        self.vmid_cgs = self.vmid.in_cgs().v
        self.dv = self.vbins[1]-self.vbins[0]
        self.dv_cgs = self.dv.in_cgs().v

        self.current_v = 0.0

        _vlos = create_vlos(normal, self.no_shifting)
        self.ds.add_field(("gas","v_los"), function=_vlos, units="cm/s")

        _intensity = self._create_intensity()
        self.ds.add_field(("gas","intensity"), function=_intensity, units=self.field_units)

        if method == "integrate" and weight_field is None:
            self.proj_units = str(ds.quan(1.0, self.field_units+"*cm").units)
        elif method == "sum":
            self.proj_units = self.field_units

        storage = {}
        pbar = get_pbar("Generating cube.", self.nv)
        for sto, i in parallel_objects(range(self.nv), storage=storage):
            self.current_v = self.vmid_cgs[i]
            if isinstance(normal, string_types):
                prj = ds.proj("intensity", ds.coordinates.axis_id[normal], method=method,
                              weight_field=weight_field, data_source=data_source)
                buf = prj.to_frb(width, self.nx, center=self.center)["intensity"]
            else:
                if data_source is None:
                    source = ds
                else:
                    source = data_source
                buf = off_axis_projection(source, self.center, normal, width,
                                          (self.nx, self.ny), "intensity",
                                          north_vector=north_vector, no_ghost=no_ghost,
                                          method=method, weight=weight_field)
            sto.result_id = i
            sto.result = buf.swapaxes(0,1)
            pbar.update(i)
        pbar.finish()

        self.data = ds.arr(np.zeros((self.nx,self.ny,self.nv)), self.proj_units)
        if is_root():
            for i, buf in sorted(storage.items()):
                self.data[:,:,i] = buf.transpose()

        self.axis_type = "velocity"

        # Now fix the width
        if iterable(self.width):
            self.width = ds.quan(self.width[0], self.width[1])
        elif not isinstance(self.width, YTQuantity):
            self.width = ds.quan(self.width, "code_length")

        self.ds.field_info.pop(("gas","intensity"))
        self.ds.field_info.pop(("gas","v_los"))

    def transform_spectral_axis(self, rest_value, units):
        """
        Change the units of the spectral axis to some equivalent unit, such
        as energy, wavelength, or frequency, by providing a *rest_value* and the
        *units* of the new spectral axis. This corresponds to the Doppler-shifting
        of lines due to gas motions and thermal broadening.
        """
        if self.axis_type != "velocity":
            self.reset_spectral_axis()
        x0 = self.ds.quan(rest_value, units)
        if x0.units.dimensions == ytdims.rate or x0.units.dimensions == ytdims.energy:
            self.vbins = x0*(1.-self.vbins.in_cgs()/clight)
        elif x0.units.dimensions == ytdims.length:
            self.vbins = x0/(1.-self.vbins.in_cgs()/clight)
        self.vmid = 0.5*(self.vbins[1:]+self.vbins[:-1])
        self.dv = self.vbins[1]-self.vbins[0]
        dims = self.dv.units.dimensions
        if dims == ytdims.rate:
            self.axis_type = "frequency"
        elif dims == ytdims.length:
            self.axis_type = "wavelength"
        elif dims == ytdims.energy:
            self.axis_type = "energy"
        elif dims == ytdims.velocity:
            self.axis_type = "velocity"

    def reset_spectral_axis(self):
        """
        Reset the spectral axis to the original velocity range and units.
        """
        self.vbins = self._vbins.copy()
        self.vmid = 0.5*(self.vbins[1:]+self.vbins[:-1])
        self.dv = self.vbins[1]-self.vbins[0]

    @parallel_root_only
    def write_fits(self, filename, clobber=False, length_unit=None,
                   sky_scale=None, sky_center=None):
        r""" Write the PPVCube to a FITS file.

        Parameters
        ----------
        filename : string
            The name of the file to write to. 
        clobber : boolean, optional
            Whether to overwrite a file with the same name that already 
            exists. Default False.
        length_unit : string, optional
            The units to convert the coordinates to in the file.
        sky_scale : tuple, optional
            Conversion between an angle unit and a length unit, if sky
            coordinates are desired, e.g. (1.0, "arcsec/kpc")
        sky_center : tuple, optional
            The (RA, Dec) coordinate in degrees of the central pixel. Must
            be specified with *sky_scale*.

        Examples
        --------
        >>> cube.write_fits("my_cube.fits", clobber=False, 
        ...                 sky_scale=(1.0,"arcsec/kpc"), sky_center=(30.,45.))
        """
        vunit = fits_info[self.axis_type][0]
        vtype = fits_info[self.axis_type][1]

        v_center = 0.5*(self.vbins[0]+self.vbins[-1]).in_units(vunit).value

        if length_unit is None:
            units = str(self.ds.get_smallest_appropriate_unit(self.width))
        else:
            units = length_unit
        units = sanitize_fits_unit(units)
        dx = self.width.in_units(units).v/self.nx
        dy = self.width.in_units(units).v/self.ny
        dv = self.dv.in_units(vunit).v

        w = _astropy.pywcs.WCS(naxis=3)
        w.wcs.crpix = [0.5*(self.nx+1), 0.5*(self.ny+1), 0.5*(self.nv+1)]
        w.wcs.cdelt = [dx,dy,dv]
        w.wcs.crval = [0.0,0.0,v_center]
        w.wcs.cunit = [units,units,vunit]
        w.wcs.ctype = ["LINEAR","LINEAR",vtype]

        fib = FITSImageData(self.data.transpose(), fields=self.field, wcs=w)
        fib.update_all_headers("bunit", re.sub('()', '', str(self.proj_units)))
        fib.update_all_headers("btype", self.field)
        if sky_scale is not None and sky_center is not None:
            fib.create_sky_wcs(sky_center, sky_scale)
        fib.writeto(filename, clobber=clobber)

    def __repr__(self):
        return "PPVCube [%d %d %d] (%s < %s < %s)" % (self.nx, self.ny, self.nv,
                                                      self.vbins[0],
                                                      fits_info[self.axis_type][2],
                                                      self.vbins[-1])

    def __getitem__(self, item):
        return self.data[item]

    def _create_intensity(self):
        def _intensity(field, data):
            v = self.current_v-data["v_los"].in_cgs().v
            T = (data["temperature"]).in_cgs().v
            w = ppv_utils.compute_weight(self.thermal_broad, self.dv_cgs,
                                         self.particle_mass, v.flatten(), T.flatten())
            w[np.isnan(w)] = 0.0
            return data[self.field]*w.reshape(v.shape)
        return _intensity

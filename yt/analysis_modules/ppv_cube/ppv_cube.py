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
from yt.utilities.fits_image import FITSImageBuffer
from yt.visualization.volume_rendering.camera import off_axis_projection
from yt.funcs import get_pbar
from yt.utilities.physical_constants import clight, mh, kboltz
import yt.units.dimensions as ytdims
from yt.units.yt_array import YTQuantity
from yt.funcs import iterable
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_root_only

def create_vlos(z_hat):
    def _v_los(field, data):
        vz = data["velocity_x"]*z_hat[0] + \
             data["velocity_y"]*z_hat[1] + \
             data["velocity_z"]*z_hat[2]
        return -vz
    return _v_los

fits_info = {"velocity":("m/s","VELOCITY","v"),
             "frequency":("Hz","FREQUENCY","f"),
             "energy":("eV","ENERGY","E"),
             "wavelength":("angstrom","WAVELENG","lambda")}

class PPVCube(object):
    def __init__(self, ds, normal, field, center="c", width=(1.0,"unitary"),
                 dims=(100,100,100), velocity_bounds=None, rest_value=None,
                 thermal_broad=False, particle_weight=56.):
        r""" Initialize a PPVCube object.

        Parameters
        ----------
        ds : dataset
            The dataset.
        normal : array_like or string
            The normal vector along with to make the projections. If an array, it
            will be normalized. If a string, it will be assumed to be along one of the
            principal axes of the domain ("x","y", or "z").
        field : string
            The field to project.
        center : float, tuple, or string

        width : float or tuple, optional
            The width of the projection in length units. Specify a float
            for code_length units or a tuple (value, units).
        dims : tuple, optional
            A 3-tuple of dimensions (nx,ny,nv) for the cube.
        velocity_bounds : tuple, optional
            A 3-tuple of (vmin, vmax, units) for the velocity bounds to
            integrate over. If None, the largest velocity of the
            dataset will be used, e.g. velocity_bounds = (-v.max(), v.max())
        rest_value : tuple or YTQuantity, optional
            A (value, unit) tuple or YTQuantity indicating the rest value of the spectral
            axis of the PPV cube. The spectral axis will be converted to the units and
            displaced by the value. Can be in units of energy, wavelength, or frequency.
            If not set the default is to leave the spectral axis in velocity units.
        particle_weight : float, optional
            Set this value to the atomic weight of the particle that is emitting the line
            if *thermal_broad* is True. Defaults to 56 (Fe).

        Examples
        --------
        >>> i = 60*np.pi/180.
        >>> L = [0.0,np.sin(i),np.cos(i)]
        >>> cube = PPVCube(ds, L, "density", width=(10.,"kpc"),
        ...                velocity_bounds=(-5.,4.,"km/s"))
        """

        self.ds = ds
        self.field = field
        self.width = width
        if isinstance(rest_value, YTQuantity) or rest_value is None:
            self.rest_value = rest_value
        else:
            self.rest_value = ds.quan(rest_value[0], rest_value[1])
        self.particle_mass = particle_weight*mh
        self.thermal_broad = thermal_broad

        self.center = ds.coordinates.sanitize_center(center, normal)[0]

        self.nx = dims[0]
        self.ny = dims[1]
        self.nv = dims[2]

        if isinstance(normal, basestring):
            los_vec = np.zeros(3)
            los_vec[ds.coordinates.axis_id[normal]] = 1.0
        else:
            orient = Orientation(normal)
            los_vec = orient.unit_vectors[2]

        dd = ds.all_data()

        fd = dd._determine_fields(field)[0]

        self.field_units = ds._get_field_info(fd).units

        if velocity_bounds is None:
            vmin, vmax = dd.quantities.extrema("velocity_magnitude")
            self.v_bnd = -vmax, vmax
        else:
            self.v_bnd = (ds.quan(velocity_bounds[0], velocity_bounds[2]),
                          ds.quan(velocity_bounds[1], velocity_bounds[2]))

        self.vbins = np.linspace(self.v_bnd[0], self.v_bnd[1], num=self.nv+1)
        self.vmid = 0.5*(self.vbins[1:]+self.vbins[:-1])
        self.vmid_cgs = self.vmid.in_cgs()
        self.dv = self.vbins[1]-self.vbins[0]
        self.dv_cgs = self.dv.in_cgs()

        _vlos = create_vlos(los_vec)
        self.ds.field_info.add_field(("gas","v_los"), function=_vlos, units="cm/s")

        if thermal_broad:
            self.v2_th = lambda T: 2.*kboltz*T/self.particle_mass
            self.phi_th = lambda v, T: self.dv_cgs*np.exp(-v*v/self.v2_th(T))/(np.sqrt(np.pi*self.v2_th(T)))
        else:
            self.v2_th = lambda T: 1.0
            self.phi_th = lambda v, T: np.maximum(1.-np.abs(v)/self.dv_cgs,0.0)

        self.data = ds.arr(np.zeros((self.nx,self.ny,self.nv)), self.field_units)
        pbar = get_pbar("Generating cube.", self.nv)
        for i in xrange(self.nv):
            _intensity = self._create_intensity(i)
            ds.add_field(("gas","intensity"), function=_intensity, units=self.field_units)
            if isinstance(normal, basestring):
                prj = ds.proj("intensity", ds.coordinates.axis_id[normal])
                buf = prj.to_frb(width, self.nx, center=self.center)["intensity"]
            else:
                buf = off_axis_projection(ds, self.center, normal, width,
                                          (self.nx, self.ny), "intensity")
            self.data[:,:,i] = buf[:,:]
            ds.field_info.pop(("gas","intensity"))
            pbar.update(i)
        pbar.finish()

        self.proj_units = self.data.units

        if self.rest_value is not None:
            # If we want units other than velocity, we re-calculate these quantities
            self.vbins = self.rest_value*(1.-self.vbins.in_cgs()/clight)
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

        # Now fix the width
        if iterable(self.width):
            self.width = ds.quan(self.width[0], self.width[1])
        else:
            self.width = ds.quan(self.width, "code_length")

    @parallel_root_only
    def write_fits(self, filename, clobber=True, sky_scale=None,
                   sky_center=None):
        r""" Write the PPVCube to a FITS file.

        Parameters
        ----------
        filename : string
            The name of the file to write.
        clobber : boolean
            Whether or not to clobber an existing file with the same name.
        sky_scale : tuple or YTQuantity
            Conversion between an angle unit and a length unit, if sky
            coordinates are desired.
            Examples: (1.0, "arcsec/kpc"), YTQuantity(0.001, "deg/kpc")
        sky_center : tuple, optional
            The (RA, Dec) coordinate in degrees of the central pixel if
            *sky_scale* has been specified.

        Examples
        --------
        >>> cube.write_fits("my_cube.fits", clobber=False, sky_scale=(1.0,"arcsec/kpc"))
        """
        if sky_scale is None:
            center = (0.0,0.0)
            types = ["LINEAR","LINEAR"]
        else:
            if iterable(sky_scale):
                sky_scale = self.ds.quan(sky_scale[0], sky_scale[1])
            if sky_center is None:
                center = (30.,45.)
            else:
                center = sky_center
            types = ["RA---TAN","DEC--TAN"]

        vunit = fits_info[self.axis_type][0]
        vtype = fits_info[self.axis_type][1]

        v_center = 0.5*(self.vbins[0]+self.vbins[-1]).in_units(vunit).value

        if sky_scale:
            dx = (self.width*sky_scale).in_units("deg")/self.nx
            units = "deg"
        else:
            dx = self.width/self.nx
            units = str(self.width.units)
        # Hack because FITS is stupid and doesn't understand case
        if units == "Mpc": units = "kpc"
        dy = dx
        dv = self.dv.in_units(vunit)

        if sky_scale:
            dx *= -1.

        w = _astropy.pywcs.WCS(naxis=3)
        w.wcs.crpix = [0.5*(self.nx+1), 0.5*(self.ny+1), 0.5*(self.nv+1)]
        w.wcs.cdelt = [dx.v,dy.v,dv.v]
        w.wcs.crval = [center[0],center[1],v_center]
        w.wcs.cunit = [units,units,vunit]
        w.wcs.ctype = [types[0],types[1],vtype]

        fib = FITSImageBuffer(self.data.transpose(), fields=self.field, wcs=w)
        fib[0].header["bunit"] = str(self.proj_units)
        fib[0].header["btype"] = self.field

        fib.writeto(filename, clobber=clobber)

    def _create_intensity(self, i):
        def _intensity(field, data):
            w = self.phi_th(self.vmid_cgs[i]-data["v_los"], data["temperature"])
            return data[self.field]*w
        return _intensity

    def __repr__(self):
        return "PPVCube [%d %d %d] (%s < %s < %s)" % (self.nx, self.ny, self.nv,
                                                      self.vbins[0],
                                                      fits_info[self.axis_type][2],
                                                      self.vbins[-1])

    def __getitem__(self, item):
        return self.data[item]

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

def create_vlos(z_hat):
    def _v_los(field, data):
        vz = data["velocity_x"]*z_hat[0] + \
             data["velocity_y"]*z_hat[1] + \
             data["velocity_z"]*z_hat[2]
        return -vz
    return _v_los

class PPVCube(object):
    def __init__(self, ds, normal, field, width=(1.0,"unitary"),
                 dims=(100,100,100), velocity_bounds=None):
        r""" Initialize a PPVCube object.

        Parameters
        ----------
        ds : dataset
            The dataset.
        normal : array_like
            The normal vector along with to make the projections.
        field : string
            The field to project.
        width : float or tuple, optional
            The width of the projection in length units. Specify a float
            for code_length units or a tuple (value, units).
        dims : tuple, optional
            A 3-tuple of dimensions (nx,ny,nv) for the cube.
        velocity_bounds : tuple, optional
            A 3-tuple of (vmin, vmax, units) for the velocity bounds to
            integrate over. If None, the largest velocity of the
            dataset will be used, e.g. velocity_bounds = (-v.max(), v.max())

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

        self.nx = dims[0]
        self.ny = dims[1]
        self.nv = dims[2]

        normal = np.array(normal)
        normal /= np.sqrt(np.dot(normal, normal))
        vecs = np.identity(3)
        t = np.cross(normal, vecs).sum(axis=1)
        ax = t.argmax()
        north = np.cross(normal, vecs[ax,:]).ravel()
        orient = Orientation(normal, north_vector=north)

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
        self.dv = (self.v_bnd[1]-self.v_bnd[0])/self.nv

        _vlos = create_vlos(orient.unit_vectors[2])
        ds.field_info.add_field(("gas","v_los"), function=_vlos, units="cm/s")

        self.data = ds.arr(np.zeros((self.nx,self.ny,self.nv)), self.field_units)
        pbar = get_pbar("Generating cube.", self.nv)
        for i in xrange(self.nv):
            _intensity = self._create_intensity(i)
            ds.add_field(("gas","intensity"), function=_intensity, units=self.field_units)
            prj = off_axis_projection(ds, ds.domain_center, normal, width,
                                      (self.nx, self.ny), "intensity")
            self.data[:,:,i] = prj[:,:]
            ds.field_info.pop(("gas","intensity"))
            pbar.update(i)

        pbar.finish()

    def write_fits(self, filename, clobber=True, length_unit=(10.0, "kpc"),
                   sky_center=(30.,45.)):
        r""" Write the PPVCube to a FITS file.

        Parameters
        ----------
        filename : string
            The name of the file to write.
        clobber : boolean
            Whether or not to clobber an existing file with the same name.
        length_unit : tuple, optional
            The length that corresponds to the width of the projection in
            (value, unit) form. Accepts a length unit or 'deg'.
        sky_center : tuple, optional
            The (RA, Dec) coordinate in degrees of the central pixel if
            *length_unit* is 'deg'.

        Examples
        --------
        >>> cube.write_fits("my_cube.fits", clobber=False, length_unit=(5,"deg"))
        """
        if length_unit[1] == "deg":
            center = sky_center
            types = ["RA---SIN","DEC--SIN"]
        else:
            center = [0.0,0.0]
            types = ["LINEAR","LINEAR"]

        v_center = 0.5*(self.v_bnd[0]+self.v_bnd[1]).in_units("m/s").value

        dx = length_unit[0]/self.nx
        dy = length_unit[0]/self.ny
        dv = self.dv.in_units("m/s").value

        if length_unit[1] == "deg":
            dx *= -1.

        w = _astropy.pywcs.WCS(naxis=3)
        w.wcs.crpix = [0.5*(self.nx+1), 0.5*(self.ny+1), 0.5*(self.nv+1)]
        w.wcs.cdelt = [dx,dy,dv]
        w.wcs.crval = [center[0], center[1], v_center]
        w.wcs.cunit = [length_unit[1],length_unit[1],"m/s"]
        w.wcs.ctype = [types[0],types[1],"VELO-LSR"]

        fib = FITSImageBuffer(self.data.transpose(), fields=self.field, wcs=w)
        fib[0].header["bunit"] = self.field_units
        fib[0].header["btype"] = self.field

        fib.writeto(filename, clobber=clobber)

    def _create_intensity(self, i):
        def _intensity(field, data):
            vlos = data["v_los"]
            w = np.abs(vlos-self.vmid[i])/self.dv.in_units(vlos.units)
            w = 1.-w
            w[w < 0.0] = 0.0
            return data[self.field]*w
        return _intensity

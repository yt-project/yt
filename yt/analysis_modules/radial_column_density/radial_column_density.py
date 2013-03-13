"""
Calculate the radial column density around a point.

Author: Stephen Skory <s@skory.us>
Affiliation: CU Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Stephen Skory.  All Rights Reserved.

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

from yt.mods import *
import yt.visualization.volume_rendering.camera as camera
import yt.utilities.lib as au
from yt.utilities.math_utils import periodic_dist
from yt.data_objects.field_info_container import FieldDetector
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface

def _col_dens(field, data):
    return data[field]

class RadialColumnDensity(ParallelAnalysisInterface):
    r"""
    Calculate radial column densities in preparation to
    adding them as a derived field.
    
    This class is the first step in calculating a derived radial 
    column density field.
    Given a central point, this calculates the column density to all cell
    centers within the given radius in units of centimeters times
    the units of the basis field.
    For example, if the basis field is `NumberDensity`, which has units
    of 1 / cm^3, the units of the derived field will be 1 / cm^2.
    Please see the documentation or the example below on how to
    use this to make the derived field which can be used like any other
    derived field.
    
    This builds a number of spherical 
    surfaces where the column density is calculated
    using HEALPix Volume Rendering. The values of the column density at
    grid points is then linearly interpolated between the two nearest
    surfaces (one inward, one outward).
    Please see the HEALPix Volume Rendering documentation for more on
    that part of this calculation.
    
    Parameters
    ----------
    pf : `StaticOutput`
        The dataset to operate on.
    field : string
        The name of the basis field over which to
        calculate a column density.
    center : array_like
        A list or array giving the location of where to
        calculate the start of
        the column density.
        This will probably be "an object of interest" like
        a star, black hole, or the center of a galaxy.
    max_radius : float
        How far out to calculate the column density, in code units. This
        value will be automatically reduced if the supplied value would
        result in calculating column densities outside the volume.
        Default = 0.5.
    steps : integer
        How many surfaces to use. A higher number is more accurate, but
        takes more resources.
        Default = 10
    base : string
        How to evenly space the surfaces: linearly with "lin" or
        logarithmically with "log".
        Default = "lin".
    Nside : int
        The resolution of column density calculation as performed by
        HEALPix. Higher numbers mean higher quality. Max = 8192.
        Default = 32.
    ang_divs : imaginary integer
        This number controls the gridding of the HEALPix projection onto
        the spherical surfaces. Higher numbers mean higher quality.
        Default = 800j.
    
    Examples
    --------
    
    >>> rcdnumdens = RadialColumnDensity(pf, 'NumberDensity', [0.5, 0.5, 0.5])
    >>> def _RCDNumberDensity(field, data, rcd = rcdnumdens):
            return rcd._build_derived_field(data)
    >>> add_field('RCDNumberDensity', _RCDNumberDensity, units=r'1/\rm{cm}^2')
    """
    def __init__(self, pf, field, center, max_radius = 0.5, steps = 10,
            base='lin', Nside = 32, ang_divs = 800j):
        ParallelAnalysisInterface.__init__(self)
        self.pf = pf
        self.center = np.asarray(center)
        self.max_radius = max_radius
        self.steps = steps
        self.base = base
        self.Nside = Nside
        self.ang_divs = ang_divs
        self.real_ang_divs = int(np.abs(ang_divs))
        self.phi, self.theta = np.mgrid[0.0:2*np.pi:ang_divs, 0:np.pi:ang_divs]
        self.phi1d = self.phi[:,0]
        self.theta1d = self.theta[0,:]
        self.dphi = self.phi1d[1] - self.phi1d[0]
        self.dtheta = self.theta1d[1] - self.theta1d[0]
        self.pixi = au.arr_ang2pix_nest(Nside, self.theta.ravel(), self.
            phi.ravel())
        self.dw = pf.domain_right_edge - pf.domain_left_edge
        # Here's where we actually do stuff.
        self._fix_max_radius()
        self._make_bins()
        self._build_surfaces(field)
    
    def _fix_max_radius(self):
        # The max_radius can only be the distance from the center point to
        # the closest face of the volume. This is because the column density
        # for a surface outside the volume is ill-defined due to the way 
        # normalization is handled in the volume render.
        # It may be possible to fix this in
        # the future, and allow these calculations in the whole volume,
        # but this will work for now.
        right = self.pf.domain_right_edge - self.center
        left = self.center - self.pf.domain_left_edge
        min_r = np.min(right)
        min_l = np.min(left)
        self.max_radius = np.min([self.max_radius, min_r, min_l])
    
    def _make_bins(self):
        # We'll make the bins start from the smallest cell size to the
        # specified radius. Column density inside the same cell as our 
        # center is kind of ill-defined, anyway.
        if self.base == 'lin':
            self.bins = np.linspace(self.pf.h.get_smallest_dx(), self.max_radius,
                self.steps)
        elif self.base == 'log':
            self.bins = np.logspace(np.log10(self.pf.h.get_smallest_dx()),
                np.log10(self.max_radius), self.steps)
    
    def _build_surfaces(self, field):
        # This will be index by bin index.
        self.surfaces = {}
        for i, radius in enumerate(self.bins):
            cam = camera.HEALpixCamera(self.center, radius, self.Nside,
                pf = self.pf, log_fields = [False], fields = field)
            bitmap = cam.snapshot()
            self.surfaces[i] = radius * self.pf['cm'] * \
                bitmap[:,0,0][self.pixi].reshape((self.real_ang_divs,self.real_ang_divs))
            self.surfaces[i] = self.comm.mpi_allreduce(self.surfaces[i], op='max')

    def _build_derived_field(self, data, minval=None):
        r"""
        Parameters
        ----------
        
        minval : float
            This parameter will set any values of the
            field that are zero to this minimum value.
            Values of zero are found outside the maximum radius and
            in the cell of the user-specified center point.
            This setting is useful if the field is going to be logged
            (e.g. np.log10) where zeros are inconvenient.
            Default = None
        """
        x = data['x']
        sh = x.shape
        ad = np.prod(sh)
        if type(data) == type(FieldDetector()):
            return np.ones(sh)
        y = data['y']
        z = data['z']
        pos = np.array([x.reshape(ad), y.reshape(ad), z.reshape(ad)]).T
        del x, y, z
        vals = self._interpolate_value(pos)
        del pos
        if minval:
            zeros = (vals == 0.)
            vals[zeros] = minval
            del zeros
        vals.shape = sh
        return vals
    
    def _interpolate_value(self, pos):
        # Given a position, find the two surfaces it's in between,
        # and the interpolate values from the surfaces to the point
        # according to the points angle.
        # 1. Find the angle from the center point to the position.
        vec = pos - self.center
        phi = np.arctan2(vec[:, 1], vec[:, 0])
        # Convert the convention from [-pi, pi) to [0, 2pi).
        sel = (phi < 0)
        phi[sel] += 2 * np.pi
        # Find the radius.
        r = np.sqrt(np.sum(vec * vec, axis = 1))
        # Keep track of the points outside of self.max_radius, which we'll
        # handle separately before we return.
        outside = (r > self.max_radius)
        theta = np.arccos(vec[:, 2] / r)
        # 2. Find the bin for this position.
        digi = np.digitize(r, self.bins)
        # Find the values on the inner and outer surfaces.
        in_val = np.zeros_like(r)
        out_val = np.zeros_like(r)
        # These two will be used for interpolation.
        in_r = np.zeros_like(r)
        out_r = np.zeros_like(r)
        for bin in np.unique(digi):
            sel = (digi == bin)
            # Special case if we're outside the largest sphere.
            if bin == len(self.bins):
                in_val[sel] = 0.
                out_val[sel] = 0.
                # Just something non-zero so we don't get divide errors later.
                in_r[sel] = .1
                out_r[sel] = .2
                continue
            # Special case if we're inside the smallest sphere.
            elif bin == 0:
                in_val[sel] = np.zeros_like(phi[sel])
                in_r[sel] = 0.
                out_val[sel] = self._interpolate_surface_value(1,
                    phi[sel], theta[sel])
                out_r[sel] = self.bins[1]
                continue
            # General case.
            else:
                in_val[sel] = self._interpolate_surface_value(bin - 1,
                    phi[sel], theta[sel])
                in_r[sel] = self.bins[bin - 1]
                out_val[sel] = self._interpolate_surface_value(bin,
                    phi[sel], theta[sel])
                out_r[sel] = self.bins[bin]
        # Interpolate using a linear fit in column density / r space.
        val = np.empty_like(r)
        # Special case for inside smallest sphere.
        sel = (digi == 0)
        val[sel] = (1. - (out_r[sel] - r[sel]) / out_r[sel]) * out_val[sel]
        np.invert(sel, sel) # In-place operation!
        val[sel] = (out_val[sel] - in_val[sel]) / (out_r[sel] - in_r[sel]) * \
            (r[sel] - in_r[sel]) + in_val[sel]
        # Fix the things to zero that should be zero.
        val[outside] = 0.
        return val
        
    def _interpolate_surface_value(self, bin, phi, theta):
        # Given a surface bin and an angle, interpolate the value on
        # that surface to the angle.
        # 1. Find the four values closest to the angle.
        phi_bin = np.digitize(phi, self.phi1d)
        theta_bin = np.digitize(theta, self.theta1d)
        val00 = self.surfaces[bin][phi_bin - 1, theta_bin - 1]
        val01 = self.surfaces[bin][phi_bin - 1, theta_bin]
        val10 = self.surfaces[bin][phi_bin, theta_bin - 1]
        val11 = self.surfaces[bin][phi_bin, theta_bin]
        # 2. Linearly interpolate the four values to the points.
        int_val0 = (val10 - val00) / self.dphi * \
            (phi - self.phi1d[phi_bin - 1]) + val00
        int_val1 = (val11 - val01) / self.dphi * \
            (phi - self.phi1d[phi_bin - 1]) + val01
        vals = (int_val1 - int_val0) / self.dtheta * \
            (theta - self.theta1d[theta_bin - 1]) + int_val0
        return vals



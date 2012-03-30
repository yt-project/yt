"""
Export/Import of volume rendered images.

Author: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Samuel Skillman.  All Rights Reserved.

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
import h5py
try: import pyfits
except: pass
import numpy as na

from yt.funcs import *

def export_rgba(image, fn, h5=True, fits=False, ):
    """
    This function accepts an *image*, of shape (N,M,4) corresponding to r,g,b,a,
    and saves to *fn*.  If *h5* is True, then it will save in hdf5 format.  If
    *fits* is True, it will save in fits format.
    """
    if h5:
        f = h5py.File('%s.h5'%fn, "w")
        f.create_dataset("R", data=image[:,:,0])
        f.create_dataset("G", data=image[:,:,1])
        f.create_dataset("B", data=image[:,:,2])
        f.create_dataset("A", data=image[:,:,3])
        f.close()
    if fits:
        try:
            hdur = pyfits.PrimaryHDU(image[:,:,0])
            hdug = pyfits.ImageHDU(image[:,:,1])
            hdub = pyfits.ImageHDU(image[:,:,2])
            hdua = pyfits.ImageHDU(image[:,:,3])
            hdulist = pyfits.HDUList([hdur,hdug,hdub,hdua])
            hdulist.writeto('%s.fits'%fn,clobber=True)
        except: print 'You do not have pyfits, install before attempting to use fits exporter'

def import_rgba(name, h5=True):
    """
    This function will read back in an HDF5 file, as saved by export_rgba, and
    return the frames to the user.  *name* is the name of the file to be read
    in.
    """
    if h5:
        f = h5py.File(name, "r")
        r = f['R'].value
        g = f['G'].value
        b = f['B'].value
        a = f['A'].value
        f.close()
    else:
        print 'No support for fits import.'
    return na.array([r,g,b,a]).swapaxes(0,2).swapaxes(0,1)

def plot_channel(image, name, cmap='gist_heat', log=True, dex=3, zero_factor=1.0e-10, 
                 label=None, label_color='w', label_size='large'):
    """
    This function will plot a single channel. *image* is an array shaped like
    (N,M), *name* is the pefix for the output filename.  *cmap* is the name of
    the colormap to apply, *log* is whether or not the channel should be
    logged.  Additionally, you may optionally specify the minimum-value cutoff
    for scaling as *dex*, which is taken with respect to the minimum value of
    the image.  *zero_factor* applies a minimum value to all zero-valued
    elements.  Optionally, *label*, *label_color* and *label_size* may be
    specified.
    """
    import matplotlib
    import pylab
    Nvec = image.shape[0]
    image[na.isnan(image)] = 0.0
    ma = image[image>0.0].max()
    image[image==0.0] = ma*zero_factor
    if log:
        mynorm = matplotlib.colors.LogNorm(ma/(10.**dex), ma)

    pylab.clf()
    pylab.gcf().set_dpi(100)
    pylab.gcf().set_size_inches((Nvec/100.0, Nvec/100.0))
    pylab.gcf().subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0, hspace=0.0)
    mycm = pylab.cm.get_cmap(cmap)
    if log:
        pylab.imshow(image,cmap=mycm, norm=mynorm, interpolation='nearest')
    else:
        pylab.imshow(image,cmap=mycm, interpolation='nearest')
    if label is not None:
        pylab.text(20, 20,label, color = label_color, size=label_size) 
    pylab.savefig("%s_%s.png" % (name,cmap))
    pylab.clf()

def plot_rgb(image, name, label=None, label_color='w', label_size='large'):
    """
    This will plot the r,g,b channels of an *image* of shape (N,M,3) or
    (N,M,4).  *name* is the prefix of the file name, which will be supplemented
    with "_rgb.png."  *label*, *label_color* and *label_size* may also be
    specified.
    """
    import pylab
    Nvec = image.shape[0]
    image[na.isnan(image)] = 0.0
    if image.shape[2] >= 4:
        image = image[:,:,:3]
    pylab.clf()
    pylab.gcf().set_dpi(100)
    pylab.gcf().set_size_inches((Nvec/100.0, Nvec/100.0))
    pylab.gcf().subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0, hspace=0.0)
    pylab.imshow(image, interpolation='nearest')
    if label is not None:
        pylab.text(20, 20, label, color = label_color, size=label_size) 
    pylab.savefig("%s_rgb.png" % name)
    pylab.clf()

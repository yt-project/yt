"""
Export/Import of volume rendered images.

Author: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Samuel Skillman.  All Rights Reserved.

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
import numpy as na
import matplotlib; from matplotlib import pylab
from yt.extensions.volume_rendering import *
from yt.funcs import *
import h5py
try: import pyfits
except: pass

'''
export_rgba:
args: 
image - rgba numpy array of dimension NxNx4 
fn - file basename for exported files

kwargs:
h5 [True]: Use hdf5 to export the rgba channels, creates fn.h5 file.
fits [False]: Export to fits format (requires pyfits), creates fn_[r,g,b,a].fits files

output: None
'''
def export_rgba(image, fn, h5=True, fits=False, ):
    if h5:
        f = h5py.File('%s.h5'%fn, "w")
        f.create_dataset("R", data=image[:,:,0])
        f.create_dataset("G", data=image[:,:,1])
        f.create_dataset("B", data=image[:,:,2])
        f.create_dataset("A", data=image[:,:,3])
        f.close()
    if fits:
        try:
            hdu = pyfits.PrimaryHDU(image[:,:,0])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('%s_r.fits'%fn,clobber=True)
            hdu = pyfits.PrimaryHDU(image[:,:,1])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('%s_g.fits'%fn,clobber=True)
            hdu = pyfits.PrimaryHDU(image[:,:,2])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('%s_b.fits'%fn,clobber=True)
            hdu = pyfits.PrimaryHDU(image[:,:,3])
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('%s_a.fits'%fn,clobber=True)
        except: print 'You do not have pyfits, install before attempting to use fits exporter'

'''
import_rgba:
args: 
name of hdf5 file containing R,G,B,A arrays

kwargs:
h5 [True]: Use hdf5 to input the rgba channels.

returns: 
a Numpy array of shape NxNx4
'''
def import_rgba(name, h5=True):
    if h5:
        f = h5py.File(name, "r")
        r = f['R'].value
        g = f['G'].value
        b = f['B'].value
        a = f['A'].value
    else:
        print 'No support for fits import.'
    return na.array([r,g,b,a]).swapaxes(0,2).swapaxes(0,1)

'''
plot_channel:
args: 
image - single channel of an image, of shape NxN
name - prefix for image output 

kwargs:
cmap ['gist_stern'] - name of colormap to use
log [True] - use a log scale
dex [3] - Scale the image between image.max()/10**dex and image.max()
zero_factor [1.0e-10] - Set all zero values to zero_factor
label [None] - Label in upper left corner of image
label_color ['w'] - Color of label
label_size ['large'] - Size of label

returns: 
None

output:
name_cmap.png
'''
def plot_channel(image, name, cmap='gist_heat', log=True, dex=3, zero_factor=1.0e-10, 
                 label=None, label_color='w', label_size='large'):
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

'''
plot_rgb:
args: 
image - rgb channels of an image, of shape NxNx3 or NxNx4.  Values must be normalized from 0.0-1.0
name - prefix for image output 

kwargs:
label [None] - Label in upper left corner of image
label_color ['w'] - Color of label
label_size ['large'] - Size of label

returns: 
None

output:
name_rgb.png
'''
def plot_rgb(image, name, label=None, label_color='w', label_size='large'):
    Nvec = image.shape[0]
    image[na.isnan(image)] = 0.0
    if image.shape[2] == 4:
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

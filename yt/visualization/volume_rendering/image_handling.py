"""
Export/Import of volume rendered images.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np

from yt.funcs import mylog

def export_rgba(image, fn, h5=True, fits=False, ):
    """
    This function accepts an *image*, of shape (N,M,4) corresponding to r,g,b,a,
    and saves to *fn*.  If *h5* is True, then it will save in hdf5 format.  If
    *fits* is True, it will save in fits format.
    """
    if (not h5 and not fits) or (h5 and fits):
        raise ValueError("Choose either HDF5 or FITS format!")
    if h5:
        f = h5py.File('%s.h5'%fn, "w")
        f.create_dataset("R", data=image[:,:,0])
        f.create_dataset("G", data=image[:,:,1])
        f.create_dataset("B", data=image[:,:,2])
        f.create_dataset("A", data=image[:,:,3])
        f.close()
    if fits:
        from yt.utilities.fits_image import FITSImageData
        data = {}
        data["r"] = image[:,:,0]
        data["g"] = image[:,:,1]
        data["b"] = image[:,:,2]
        data["a"] = image[:,:,3]
        nx, ny = data["r"].shape
        fib = FITSImageData(data)
        fib.writeto('%s.fits'%fn,clobber=True)

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
        mylog.error('No support for fits import.')
    return np.array([r,g,b,a]).swapaxes(0,2).swapaxes(0,1)

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
    image[np.isnan(image)] = 0.0
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
    image[np.isnan(image)] = 0.0
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

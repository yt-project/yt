"""
ImageArray Class

Authors: Samuel Skillman <samskillman@gmail.com>
Affiliation: University of Colorado at Boulder

Homepage: http://yt-project.org/
License:
    Copyright (C) 2012 Samuel Skillman.  All Rights Reserved.

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

import numpy as np
import h5py as h5
from yt.visualization.image_writer import write_bitmap, write_image

class ImageArray(np.ndarray):
    r"""A custom Numpy ndarray used for images.

    This differs from ndarray in that you can optionally specify an
    info dictionary which is used later in saving, and can be accessed with
    ImageArray.info.

    Parameters
    ----------
    input_array: array_like
        A numpy ndarray, or list. 

    Other Parameters
    ----------------
    info: dictionary
        Contains information to be stored with image.

    Returns
    -------
    obj: ImageArray object 

    Raises
    ------
    None

    See Also
    --------
    numpy.ndarray : Inherits

    Notes
    -----

    References
    ----------

    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.  Use the variables 'pf' for the parameter file, 'pc' for
    a plot collection, 'c' for a center, and 'L' for a vector. 

    >>> im = np.zeros([64,128,3])
    >>> for i in xrange(im.shape[0]):
    >>>     for k in xrange(im.shape[2]):
    >>>         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

    >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
    >>>     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
    >>>     'width':0.245, 'units':'cm', 'type':'rendering'}

    >>> im_arr = ImageArray(im, info=myinfo)
    >>> im_arr.save('test_ImageArray')

    Numpy ndarray documentation appended:

    """
    def __new__(cls, input_array, info=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        if info is None:
            info = {}
        obj.info = info
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    def write_hdf5(self, filename):
        r"""Writes ImageArray to hdf5 file.

        Parameters
        ----------
        filename: string
            Note filename not be modified.
       
        Examples
        -------- 
        >>> im = np.zeros([64,128,3])
        >>> for i in xrange(im.shape[0]):
        >>>     for k in xrange(im.shape[2]):
        >>>         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

        >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        >>>     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        >>>     'width':0.245, 'units':'cm', 'type':'rendering'}

        >>> im_arr = ImageArray(im, info=myinfo)
        >>> im_arr.write_hdf5('test_ImageArray.h5')

        """
        array_name = self.info.get("name","image")

        f = h5.File(filename)
        if array_name in f.keys():
            del f[array_name]
        d = f.create_dataset(array_name, data=self)
        for k, v in self.info.iteritems():
            d.attrs.create(k, v)
        f.close()

    def write_png(self, filename, clip_ratio=None):
        r"""Writes ImageArray to png file.

        Parameters
        ----------
        filename: string
            Note filename not be modified.
       
        Examples
        --------
        
        >>> im = np.zeros([64,128,3])
        >>> for i in xrange(im.shape[0]):
        >>>     for k in xrange(im.shape[2]):
        >>>         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

        >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        >>>     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        >>>     'width':0.245, 'units':'cm', 'type':'rendering'}

        >>> im_arr = ImageArray(im, info=myinfo)
        >>> im_arr.write_png('test_ImageArray.png')

        """
        if filename[-4:] != '.png': 
            filename += '.png'

        if clip_ratio is not None:
            return write_bitmap(self.swapaxes(0, 1), filename,
                                clip_ratio * self.std())
        else:
            return write_bitmap(self.swapaxes(0, 1), filename)

    def write_image(self, filename, color_bounds=None, channel=None,  cmap_name="algae", func=lambda x: x):
        r"""Writes a single channel of the ImageArray to a png file.

        Parameters
        ----------
        filename: string
            Note filename not be modified.
       
        Other Parameters
        ----------------
        channel: int
            Which channel to write out as an image. Defaults to 0
        cmap_name: string
            Name of the colormap to be used.
        color_bounds : tuple of floats, optional
            The min and max to scale between.  Outlying values will be clipped.
        cmap_name : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        func : function, optional
            A function to transform the buffer before applying a colormap. 

        Returns
        -------
        scaled_image : uint8 image that has been saved
        
        Examples
        --------
        
        >>> im = np.zeros([64,128])
        >>> for i in xrange(im.shape[0]):
        >>>     im[i,:] = np.linspace(0.,0.3*k, im.shape[1])

        >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        >>>     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        >>>     'width':0.245, 'units':'cm', 'type':'rendering'}

        >>> im_arr = ImageArray(im, info=myinfo)
        >>> im_arr.write_image('test_ImageArray.png')

        """
        if filename[-4:] != '.png': 
            filename += '.png'

        if channel is None:
            return write_image(self.swapaxes(0,1), filename, 
                               color_bounds=color_bounds, cmap_name=cmap_name, 
                               func=func)
        else:
            return write_image(self.swapaxes(0,1)[:,:,channel], filename, 
                               color_bounds=color_bounds, cmap_name=cmap_name, 
                               func=func)

    def save(self, filename, png=True, hdf5=True):
        """
        Saves ImageArray. 

        Arguments:
          filename: string
            This should not contain the extension type (.png, .h5, ...)

        Optional Arguments:
          png: boolean, default True
            Save to a png

          hdf5: boolean, default True
            Save to hdf5 file, including info dictionary as attributes.

        """
        if png:
            if len(self.shape) > 2:
                self.write_png("%s.png" % filename)
            else:
                self.write_image("%s.png" % filename)
        if hdf5:
            self.write_hdf5("%s.h5" % filename)

    __doc__ += np.ndarray.__doc__

if __name__ == "__main__":
    im = np.zeros([64,128,3])
    for i in xrange(im.shape[0]):
        for k in xrange(im.shape[2]):
            im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

    myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        'width':0.245, 'units':'cm', 'type':'rendering'}

    im_arr = ImageArray(im, info=myinfo)
    im_arr.save('test_3d_ImageArray')

    im = np.zeros([64,128])
    for i in xrange(im.shape[0]):
        im[i,:] = np.linspace(0.,0.3*k, im.shape[1])

    myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]), 
        'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),  
        'width':0.245, 'units':'cm', 'type':'rendering'}

    im_arr = ImageArray(im, info=myinfo)
    im_arr.save('test_2d_ImageArray')


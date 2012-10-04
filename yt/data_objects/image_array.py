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
from yt.visualization.image_writer import write_bitmap

class ImageArray(np.ndarray):
    """
    A custom Numpy ndarray used for images.

    This differs from ndarray in that you can optionally specify an
    info dictionary which is used later in saving, and can be accessed with
    ImageArray.info.

    Optional Arguments:
        info: dictionary
        Contains information to be stored with image.

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
        """
        Writes ImageArray to hdf5 file.

        Arguments:
            filename: string
            Note filename not be modified.

        Returns:
            None

        """
        try:
            array_name = self.info['name']
        except KeyError:
            array_name = 'image'

        f = h5.File(filename)
        if array_name in f.keys():
            del f[array_name]
        d = f.create_dataset(array_name, data=self)
        for k, v in self.info.iteritems():
            d.attrs.create(k, v)
        f.close()

    def write_png(self, filename, clip_ratio=None):
        """
        Writes ImageArray to png.

        Arguments:
            filename: string
            '.png' will be appended if not present.

        Returns:
            The bitmap array written

        Note: when writing to png, we invert the y axis
        such to prepare for the write_bitmap call.  This puts the (0,0) pixel
        in the lower left

        """
        if filename[-4:] != '.png': 
            filename += '.png'

        if clip_ratio is not None:
            return write_bitmap(self.swapaxes(0, 1), filename,
                                clip_ratio * self.std())
        else:
            return write_bitmap(self.swapaxes(0, 1), filename)

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
            self.write_png("%s.png" % filename)
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
    im_arr.save('test_ImageArray')


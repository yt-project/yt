"""
ImageArray Class

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import warnings
import numpy as np
from yt.config import \
    ytcfg
from yt.visualization.image_writer import write_bitmap, write_image
from yt.units.yt_array import YTArray


class ImageArray(YTArray):
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
    use the function.  Use the variables 'ds' for the dataset, 'pc' for
    a plot collection, 'c' for a center, and 'L' for a vector.

    >>> im = np.zeros([64,128,3])
    >>> for i in range(im.shape[0]):
    ...     for k in range(im.shape[2]):
    ...         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

    >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]),
    ...     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),
    ...     'width':0.245, 'units':'cm', 'type':'rendering'}

    >>> im_arr = ImageArray(im, info=myinfo)
    >>> im_arr.save('test_ImageArray')

    Numpy ndarray documentation appended:

    """
    def __new__(cls, input_array, input_units=None, registry=None, info=None):
        obj = super(ImageArray, cls).__new__(cls, input_array, input_units, registry)
        if info is None:
            info = {}
        obj.info = info
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        super(ImageArray, self).__array_finalize__(obj)
        self.info = getattr(obj, 'info', None)

    def write_hdf5(self, filename, dataset_name=None):
        r"""Writes ImageArray to hdf5 file.

        Parameters
        ----------
        filename: string
        The filename to create and write a dataset to

        dataset_name: string
            The name of the dataset to create in the file.

        Examples
        --------
        >>> im = np.zeros([64,128,3])
        >>> for i in range(im.shape[0]):
        ...     for k in range(im.shape[2]):
        ...         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

        >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]),
        ...     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),
        ...     'width':0.245, 'units':'cm', 'type':'rendering'}

        >>> im_arr = ImageArray(im, info=myinfo)
        >>> im_arr.write_hdf5('test_ImageArray.h5')

        """
        if dataset_name is None:
            dataset_name = self.info.get("name", "image")
        super(ImageArray, self).write_hdf5(filename, dataset_name=dataset_name,
                                           info=self.info)

    def add_background_color(self, background='black', inline=True):
        r"""Adds a background color to a 4-channel ImageArray

        This adds a background color to a 4-channel ImageArray, by default
        doing so inline.  The ImageArray must already be normalized to the
        [0,1] range.

        Parameters
        ----------
        background:
            This can be used to set a background color for the image, and can
            take several types of values:

               * ``white``: white background, opaque
               * ``black``: black background, opaque
               * ``None``: transparent background
               * 4-element array [r,g,b,a]: arbitrary rgba setting.

            Default: 'black'
        inline: boolean, optional
            If True, original ImageArray is modified. If False, a copy is first
            created, then modified. Default: True

        Returns
        -------
        out: ImageArray
            The modified ImageArray with a background color added.

        Examples
        --------
        >>> im = np.zeros([64,128,4])
        >>> for i in range(im.shape[0]):
        ...     for k in range(im.shape[2]):
        ...         im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

        >>> im_arr = ImageArray(im)
        >>> im_arr.rescale()
        >>> new_im = im_arr.add_background_color([1.,0.,0.,1.], inline=False)
        >>> new_im.write_png('red_bg.png')
        >>> im_arr.add_background_color('black')
        >>> im_arr.write_png('black_bg.png')
        """
        assert(self.shape[-1] == 4)

        if background is None:
            background = (0., 0., 0., 0.)
        elif background == 'white':
            background = (1., 1., 1., 1.)
        elif background == 'black':
            background = (0., 0., 0., 1.)

        # Alpha blending to background
        if inline:
            out = self
        else:
            out = self.copy()

        for i in range(3):
            out[:, :, i] = self[:, :, i]*self[:, :, 3]
            out[:, :, i] += background[i]*background[3]*(1.0-self[:, :, 3])
        out[:, :, 3] = self[:, :, 3]+background[3]*(1.0-self[:, :, 3])
        return out

    def rescale(self, cmax=None, amax=None, inline=True):
        r"""Rescales the image to be in [0,1] range.

        Parameters
        ----------
        cmax: float, optional
            Normalization value to use for rgb channels. Defaults to None,
            corresponding to using the maximum value in the rgb channels.
        amax: float, optional
            Normalization value to use for alpha channel. Defaults to None,
            corresponding to using the maximum value in the alpha channel.
        inline: boolean, optional
            Specifies whether or not the rescaling is done inline. If false,
            a new copy of the ImageArray will be created, returned.
            Default:True.

        Returns
        -------
        out: ImageArray
            The rescaled ImageArray, clipped to the [0,1] range.

        Notes
        -----
        This requires that the shape of the ImageArray to have a length of 3,
        and for the third dimension to be >= 3.  If the third dimension has
        a shape of 4, the alpha channel will also be rescaled.

        Examples
        --------
        >>> im = np.zeros([64,128,4])
        >>> for i in range(im.shape[0]):
        ...     for k in range(im.shape[2]):
        ...         im[i,:,k] = np.linspace(0.,0.3*k, im.shape[1])

        >>> im = ImageArray(im)
        >>> im.write_png('original.png')
        >>> im.rescale()
        >>> im.write_png('normalized.png')

        """
        assert(len(self.shape) == 3)
        assert(self.shape[2] >= 3)
        if inline:
            out = self
        else:
            out = self.copy()
        if cmax is None:
            cmax = self[:, :, :3].sum(axis=2).max()

        np.multiply(self[:, :, :3], 1.0/cmax, out[:, :, :3])

        if self.shape[2] == 4:
            if amax is None:
                amax = self[:, :, 3].max()
            if amax > 0.0:
                np.multiply(self[:, :, 3], 1.0/amax, out[:, :, 3])

        np.clip(out, 0.0, 1.0, out)
        return out

    def write_png(self, filename, sigma_clip=None, background='black',
                  rescale=True, clip_ratio=None):
        r"""Writes ImageArray to png file.

        Parameters
        ----------
        filename: string
            Filename to save to.  If None, PNG contents will be returned as a
            string.
        sigma_clip: float, optional
            Image will be clipped before saving to the standard deviation
            of the image multiplied by this value.  Useful for enhancing
            images. Default: None
        background:
            This can be used to set a background color for the image, and can
            take several types of values:

               * ``white``: white background, opaque
               * ``black``: black background, opaque
               * ``None``: transparent background
               * 4-element array [r,g,b,a]: arbitrary rgba setting.

            Default: 'black'
        rescale: boolean, optional
            If True, will write out a rescaled image (without modifying the
            original image). Default: True

        Examples
        --------
        >>> im = np.zeros([64,128,4])
        >>> for i in range(im.shape[0]):
        ...     for k in range(im.shape[2]):
        ...         im[i,:,k] = np.linspace(0.,10.*k, im.shape[1])

        >>> im_arr = ImageArray(im)
        >>> im_arr.write_png('standard.png')
        >>> im_arr.write_png('non-scaled.png', rescale=False)
        >>> im_arr.write_png('black_bg.png', background='black')
        >>> im_arr.write_png('white_bg.png', background='white')
        >>> im_arr.write_png('green_bg.png', background=[0,1,0,1])
        >>> im_arr.write_png('transparent_bg.png', background=None)

        """
        if rescale:
            scaled = self.rescale(inline=False)
        else:
            scaled = self

        if self.shape[-1] == 4:
            out = scaled.add_background_color(background, inline=False)
        else:
            out = scaled

        if filename is not None and filename[-4:] != '.png':
            filename += '.png'

        if clip_ratio is not None:
            warnings.warn("'clip_ratio' keyword is deprecated. Use 'sigma_clip' instead")
            sigma_clip = clip_ratio

        if sigma_clip is not None:
            nz = out[:, :, :3][out[:, :, :3].nonzero()]
            return write_bitmap(out.swapaxes(0, 1), filename,
                                nz.mean() + sigma_clip * nz.std())
        else:
            return write_bitmap(out.swapaxes(0, 1), filename)

    def write_image(self, filename, color_bounds=None, channel=None,
                    cmap_name=None, func=lambda x: x):
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
        >>> for i in range(im.shape[0]):
        ...     im[i,:] = np.linspace(0.,0.3*i, im.shape[1])

        >>> myinfo = {'field':'dinosaurs', 'east_vector':np.array([1.,0.,0.]),
        ...     'north_vector':np.array([0.,0.,1.]), 'normal_vector':np.array([0.,1.,0.]),
        ...     'width':0.245, 'units':'cm', 'type':'rendering'}

        >>> im_arr = ImageArray(im, info=myinfo)
        >>> im_arr.write_image('test_ImageArray.png')

        """
        if cmap_name is None:
            cmap_name = ytcfg.get("yt", "default_colormap")
        if filename[-4:] != '.png':
            filename += '.png'

        #TODO: Write info dict as png metadata
        if channel is None:
            return write_image(self.swapaxes(0, 1).to_ndarray(), filename,
                               color_bounds=color_bounds, cmap_name=cmap_name,
                               func=func)
        else:
            return write_image(self.swapaxes(0, 1)[:, :, channel].to_ndarray(),
                               filename,
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

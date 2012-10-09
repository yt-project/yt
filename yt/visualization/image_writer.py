"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation:  UCSD
License:
  Copyright (C) 2010-2011 Matthew Turk  All Rights Reserved.

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

import types
import imp
import os
import numpy as np

from yt.funcs import *
import _colormap_data as cmd
import yt.utilities.lib as au

def scale_image(image, mi=None, ma=None):
    r"""Scale an image ([NxNxM] where M = 1-4) to be uint8 and values scaled 
    from [0,255].

    Parameters
    ----------
    image : array_like or tuple of image info

    Examples
    --------

        >>> image = scale_image(image)

        >>> image = scale_image(image, min=0, max=1000)
    """
    if isinstance(image, np.ndarray) and image.dtype == np.uint8:
        return image
    if isinstance(image, (types.TupleType, types.ListType)):
        image, mi, ma = image
    if mi is None:
        mi = image.min()
    if ma is None:
        ma = image.max()
    image = (np.clip((image-mi)/(ma-mi) * 255, 0, 255)).astype('uint8')
    return image

def multi_image_composite(fn, red_channel, blue_channel,
                          green_channel = None, alpha_channel = None):
    r"""Write an image with different color channels corresponding to different
    quantities.

    Accepts at least a red and a blue array, of shape (N,N) each, that are
    optionally scaled and composited into a final image, written into `fn`.
    Can also accept green and alpha.

    Parameters
    ----------
    fn : string
        Filename to save
    red_channel : array_like or tuple of image info
        Array, of shape (N,N), to be written into the red channel of the output
        image.  If not already uint8, will be converted (and scaled) into
        uint8.  Optionally, you can also specify a tuple that includes scaling
        information, in the form of (array_to_plot, min_value_to_scale,
        max_value_to_scale).
    blue_channel : array_like or tuple of image info
        Array, of shape (N,N), to be written into the blue channel of the output
        image.  If not already uint8, will be converted (and scaled) into
        uint8.  Optionally, you can also specify a tuple that includes scaling
        information, in the form of (array_to_plot, min_value_to_scale,
        max_value_to_scale).
    green_channel : array_like or tuple of image info, optional
        Array, of shape (N,N), to be written into the green channel of the
        output image.  If not already uint8, will be converted (and scaled)
        into uint8.  If not supplied, will be left empty.  Optionally, you can
        also specify a tuple that includes scaling information, in the form of
        (array_to_plot, min_value_to_scale, max_value_to_scale).

    alpha_channel : array_like or tuple of image info, optional
        Array, of shape (N,N), to be written into the alpha channel of the output
        image.  If not already uint8, will be converted (and scaled) into uint8.
        If not supplied, will be made fully opaque.  Optionally, you can also
        specify a tuple that includes scaling information, in the form of
        (array_to_plot, min_value_to_scale, max_value_to_scale).

    Examples
    --------

        >>> red_channel = np.log10(frb["Temperature"])
        >>> blue_channel = np.log10(frb["Density"])
        >>> multi_image_composite("multi_channel1.png", red_channel, blue_channel)

    """
    red_channel = scale_image(red_channel)
    blue_channel = scale_image(blue_channel)
    if green_channel is None:
        green_channel = np.zeros(red_channel.shape, dtype='uint8')
    else:
        green_channel = scale_image(green_channel)
    if alpha_channel is None:
        alpha_channel = np.zeros(red_channel.shape, dtype='uint8') + 255
    else:
        alpha_channel = scale_image(alpha_channel) 
    image = np.array([red_channel, green_channel, blue_channel, alpha_channel])
    image = image.transpose().copy() # Have to make sure it's contiguous 
    au.write_png(image, fn)

def write_bitmap(bitmap_array, filename, max_val = None, transpose=False):
    r"""Write out a bitmapped image directly to a PNG file.

    This accepts a three- or four-channel `bitmap_array`.  If the image is not
    already uint8, it will be scaled and converted.  If it is not four channel, a
    fourth alpha channel will be added and set to fully opaque.  The resultant
    image will be directly written to `filename` as a PNG with no colormap
    applied.  `max_val` is a value used if the array is passed in as anything
    other than uint8; it will be the value used for scaling and clipping when the
    array is converted.  Additionally, the minimum is assumed to be zero; this
    makes it primarily suited for the results of volume rendered images, rather
    than misaligned projections.

    Parameters
    ----------
    bitmap_array : array_like
        Array of shape (N,M,3) or (N,M,4), to be written.  If it is not already
        a uint8 array, it will be scaled and converted to uint8.
    filename : string
        Filename to save to.  If None, PNG contents will be returned as a
        string.
    max_val : float, optional
        The upper limit to clip values to in the output, if converting to uint8.
        If `bitmap_array` is already uint8, this will be ignore.
    """
    if bitmap_array.dtype != np.uint8:
        if max_val is None: max_val = bitmap_array.max()
        bitmap_array = np.clip(bitmap_array / max_val, 0.0, 1.0) * 255
        bitmap_array = bitmap_array.astype("uint8")
    if len(bitmap_array.shape) != 3 or bitmap_array.shape[-1] not in (3,4):
        raise RuntimeError
    if bitmap_array.shape[-1] == 3:
        s1, s2 = bitmap_array.shape[:2]
        alpha_channel = 255*np.ones((s1,s2,1), dtype='uint8')
        bitmap_array = np.concatenate([bitmap_array, alpha_channel], axis=-1)
    if transpose:
        bitmap_array = bitmap_array.swapaxes(0,1)
    if filename is not None:
        au.write_png(bitmap_array.copy(), filename)
    else:
        return au.write_png_to_string(bitmap_array.copy())
    return bitmap_array

def write_image(image, filename, color_bounds = None, cmap_name = "algae", func = lambda x: x):
    r"""Write out a floating point array directly to a PNG file, scaling it and
    applying a colormap.

    This function will scale an image and directly call libpng to write out a
    colormapped version of that image.  It is designed for rapid-fire saving of
    image buffers generated using `yt.visualization.api.FixedResolutionBuffers` and the like.

    Parameters
    ----------
    image : array_like
        This is an (unscaled) array of floating point values, shape (N,N,) to
        save in a PNG file.
    filename : string
        Filename to save as.
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

    >>> sl = pf.h.slice(0, 0.5, "Density")
    >>> frb1 = FixedResolutionBuffer(sl, (0.2, 0.3, 0.4, 0.5),
                    (1024, 1024))
    >>> write_image(frb1["Density"], "saved.png")
    """
    if len(image.shape) == 3:
        mylog.info("Using only channel 1 of supplied image")
        image = image[:,:,0]
    to_plot = apply_colormap(image, color_bounds = color_bounds, cmap_name = cmap_name)
    au.write_png(to_plot, filename)
    return to_plot

def apply_colormap(image, color_bounds = None, cmap_name = 'algae', func=lambda x: x):
    r"""Apply a colormap to a floating point image, scaling to uint8.

    This function will scale an image and directly call libpng to write out a
    colormapped version of that image.  It is designed for rapid-fire saving of
    image buffers generated using `yt.visualization.api.FixedResolutionBuffers` and the like.

    Parameters
    ----------
    image : array_like
        This is an (unscaled) array of floating point values, shape (N,N,) to
        save in a PNG file.
    color_bounds : tuple of floats, optional
        The min and max to scale between.  Outlying values will be clipped.
    cmap_name : string, optional
        An acceptable colormap.  See either yt.visualization.color_maps or
        http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
    func : function, optional
        A function to transform the buffer before applying a colormap. 

    Returns
    -------
    to_plot : uint8 image with colorbar applied.

    """
    image = func(image)
    if color_bounds is None:
        mi = np.nanmin(image[~np.isinf(image)])
        ma = np.nanmax(image[~np.isinf(image)])
        color_bounds = mi, ma
    else:
        color_bounds = [func(c) for c in color_bounds]
    image = (image - color_bounds[0])/(color_bounds[1] - color_bounds[0])
    to_plot = map_to_colors(image, cmap_name)
    to_plot = np.clip(to_plot, 0, 255)
    return to_plot

def annotate_image(image, text, xpos, ypos, font_name = "Vera",
                   font_size = 24, dpi = 100):
    r"""Add text on to an existing uint8 bitmap array.

    This function accepts an image array and then directly calls freetype to
    add text on top of that array.  No array is returned.

    Parameters
    ----------
    image : array_like
        This is a (scaled) array of UINT8 values, shape (N,N,[3,4]) to
        overplot text on.
    text : string
        Text to place
    xpos : int
        The starting point, in pixels, of the text along the x axis.
    ypos : int
        The starting point, in pixels, of the text along the y axis.  Note that
        0 will be the top of the image, not the bottom.
    font_name : string (optional)
        The font to load.
    font_size : int (optional)
        Font size in points of the overlaid text.
    dpi : int (optional)
        Dots per inch for calculating the font size in pixels.
        
    Returns
    -------
    Nothing

    Examples
    --------

    >>> sl = pf.h.slice(0, 0.5, "Density")
    >>> frb1 = FixedResolutionBuffer(sl, (0.2, 0.3, 0.4, 0.5),
                    (1024, 1024))
    >>> bitmap = write_image(frb1["Density"], "saved.png")
    >>> annotate_image(bitmap, "Hello!", 0, 100)
    >>> write_bitmap(bitmap, "saved.png")
    """
    if len(image.shape) != 3 or image.dtype != np.uint8:
        raise RuntimeError("This routine requires a UINT8 bitmapped image.")
    font_path = os.path.join(imp.find_module("matplotlib")[1],
                             "mpl-data/fonts/ttf/",
                             "%s.ttf" % font_name)
    if not os.path.isfile(font_path):
        mylog.error("Could not locate %s", font_path)
        raise IOError(font_path)
    # The hard-coded 0 is the font face index.
    au.simple_writing(font_path, 0, dpi, font_size, text, image, xpos, ypos)

def map_to_colors(buff, cmap_name):
    if cmap_name not in cmd.color_map_luts:
        print "Your color map was not found in the extracted colormap file."
        raise KeyError(cmap_name)
    lut = cmd.color_map_luts[cmap_name]
    x = np.mgrid[0.0:1.0:lut[0].shape[0]*1j]
    shape = buff.shape
    mapped = np.dstack(
            [(np.interp(buff, x, v)*255) for v in lut ]).astype("uint8")
    return mapped.copy("C")

def strip_colormap_data(fn = "color_map_data.py",
            cmaps = ("jet", "algae", "hot", "gist_stern", "RdBu",
                     "kamae")):
    import pprint
    import color_maps as rcm
    f = open(fn, "w")
    f.write("### Auto-generated colormap tables, taken from Matplotlib ###\n\n")
    f.write("from numpy import array\n")
    f.write("color_map_luts = {}\n\n\n")
    if cmaps is None: cmaps = rcm.ColorMaps
    for cmap_name in sorted(cmaps):
        print "Stripping", cmap_name
        vals = rcm._extract_lookup_table(cmap_name)
        f.write("### %s ###\n\n" % (cmap_name))
        f.write("color_map_luts['%s'] = \\\n" % (cmap_name))
        f.write("   (\n")
        for v in vals:
            f.write(pprint.pformat(v, indent=3))
            f.write(",\n")
        f.write("   )\n\n")
    f.close()

def splat_points(image, points_x, points_y,
                 contribution = None, transposed = False):
    if contribution is None:
        contribution = 100.0
    val = contribution * 1.0/points_x.size
    if transposed:
        points_y = 1.0 - points_y
        points_x = 1.0 - points_x
    im = image.copy()
    au.add_points_to_image(im, points_x, points_y, val)
    return im

def write_projection(data, filename, colorbar=True, colorbar_label=None, 
                     title=None, limits=None, take_log=True, var_fig_size=False,
                     cmap='algae'):
    r"""Write a projection or volume rendering to disk with a variety of 
    pretty parameters such as limits, title, colorbar, etc.  write_projection
    uses the standard matplotlib interface to create the figure.  N.B. This code
    only works *after* you have created the projection using the standard 
    framework (i.e. the Camera interface or off_axis_projection).

    Accepts an NxM sized array representing the projection itself as well
    as the filename to which you will save this figure.  

    Parameters
    ----------
    data : array_like 
        image array as output by off_axis_projection or camera.snapshot()
    filename : string 
        the filename where the data will be saved
    colorbar : boolean
        do you want a colorbar generated to the right of the image?
    colorbar_label : string
        the label associated with your colorbar
    title : string
        the label at the top of the figure
    limits : 2-element array_like
        the lower limit and the upper limit to be plotted in the figure 
        of the data array
    take_log : boolean
        plot the log of the data array (and take the log of the limits if set)?
    var_fig_size : boolean
        If we want the resolution (and size) of the output image to scale 
        with the resolution of the image array.  
    cmap : string
        The name of the colormap.

    Examples
    --------

    >>> image = off_axis_projection(pf, c, L, W, N, "Density", no_ghost=False)
    >>> write_projection(image, 'test.png', 
                         colorbar_label="Column Density (cm$^{-2}$)", 
                         title="Offaxis Projection", limits=(1e-3,1e-5), 
                         take_log=True)
    """
    import matplotlib
    from ._mpl_imports import FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS

    # If this is rendered as log, then apply now.
    if take_log:
        norm = matplotlib.colors.LogNorm()
    else:
        norm = matplotlib.colors.Normalize()
    
    if limits is None:
        limits = [None, None]

    # Create the figure and paint the data on
    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot(111)
    fig.tight_layout()

    cax = ax.imshow(data, vmin=limits[0], vmax=limits[1], norm=norm, cmap=cmap)
    
    if title:
        ax.set_title(title)

    # Suppress the x and y pixel counts
    ax.set_xticks(())
    ax.set_yticks(())

    # Add a color bar and label if requested
    if colorbar:
        cbar = fig.colorbar(cax)
        if colorbar_label:
            cbar.ax.set_ylabel(colorbar_label)

    # If we want the resolution of the image to scale with the resolution
    # of the image array. we increase the dpi value accordingly
    
    if var_fig_size:
        N = data.shape[0]
        mag_factor = N/480.
        dpi = 100*mag_factor
    else:
        dpi = None

    if filename[-4:] in ['.png','.pdf','.eps']:
        suffix = filename[-4:]
    elif filename[-3:] == ['.ps']:
        suffix = filename[-3:]
    else:
        suffix = '.png'
        filename = "%s%s" % (filename, suffix)
    mylog.info("Saving plot %s", filename)
    if suffix == ".png":
        canvas = FigureCanvasAgg(fig)
    elif suffix == ".pdf":
        canvas = FigureCanvasPdf(fig)
    elif suffix in (".eps", ".ps"):
        canvas = FigureCanvasPS
    canvas.print_figure(filename)
    return filename


def write_fits(image, filename_prefix, clobber=True, coords=None, gzip_file=False) :

    """
    This will export a FITS image of a floating point array. The output filename is
    *filename_prefix*. If clobber is set to True, this will overwrite any existing
    FITS file.
    
    This requires the *pyfits* module, which is a standalone module
    provided by STSci to interface with FITS-format files.
    """
    r"""Write out a floating point array directly to a FITS file, optionally
    adding coordinates. 
        
    Parameters
    ----------
    image : array_like
        This is an (unscaled) array of floating point values, shape (N,N,) to save
        in a FITS file.
    filename_prefix : string
        This prefix will be prepended to every FITS file name.
    clobber : boolean
        If the file exists, this governs whether we will overwrite.
    coords : dictionary, optional
        A set of header keys and values to write to the FITS header to set up
        a coordinate system. 
    gzip_file : boolean, optional
        gzip the file after writing, default False
    """
    
    import pyfits
    from os import system
    
    if filename_prefix.endswith('.fits'): filename_prefix=filename_prefix[:-5]
    
    hdu = pyfits.PrimaryHDU(image)

    if (coords is not None) :

        hdu.header.update('WCSNAMEP', "PHYSICAL")
        hdu.header.update('CTYPE1P', "LINEAR")
        hdu.header.update('CTYPE2P', "LINEAR")
        hdu.header.update('CRPIX1P', 0.5)
        hdu.header.update('CRPIX2P', 0.5)
        hdu.header.update('CRVAL1P', coords["xmin"])
        hdu.header.update('CRVAL2P', coords["ymin"])
        hdu.header.update('CDELT1P', coords["dx"])
        hdu.header.update('CDELT2P', coords["dy"])
        
        hdu.header.update('CTYPE1', "LINEAR")
        hdu.header.update('CTYPE2', "LINEAR")
        hdu.header.update('CUNIT1', coords["units"])
        hdu.header.update('CUNIT2', coords["units"])
        hdu.header.update('CRPIX1', 0.5)
        hdu.header.update('CRPIX2', 0.5)
        hdu.header.update('CRVAL1', coords["xmin"])
        hdu.header.update('CRVAL2', coords["ymin"])
        hdu.header.update('CDELT1', coords["dx"])
        hdu.header.update('CDELT2', coords["dy"])

    hdu.writeto("%s.fits" % (filename_prefix), clobber=clobber)

    if (gzip_file) :
        clob = ""
        if (clobber) : clob="-f"
        system("gzip "+clob+" %s.fits" % (filename_prefix))
    

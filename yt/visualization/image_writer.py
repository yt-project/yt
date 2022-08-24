import builtins

import numpy as np

from yt.config import ytcfg
from yt.funcs import mylog
from yt.units.yt_array import YTQuantity
from yt.utilities import png_writer as pw
from yt.utilities.exceptions import YTNotInsideNotebook
from yt.utilities.lib import image_utilities as au
from yt.visualization.color_maps import get_colormap_lut

from ._commons import get_canvas, validate_image_name


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
    if isinstance(image, (tuple, list)):
        image, mi, ma = image
    if mi is None:
        mi = image.min()
    if ma is None:
        ma = image.max()
    image = (np.clip((image - mi) / (ma - mi) * 255, 0, 255)).astype("uint8")
    return image


def multi_image_composite(
    fn, red_channel, blue_channel, green_channel=None, alpha_channel=None
):
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

        >>> red_channel = np.log10(frb[("gas", "temperature")])
        >>> blue_channel = np.log10(frb[("gas", "density")])
        >>> multi_image_composite("multi_channel1.png", red_channel, blue_channel)

    """
    red_channel = scale_image(red_channel)
    blue_channel = scale_image(blue_channel)
    if green_channel is None:
        green_channel = np.zeros(red_channel.shape, dtype="uint8")
    else:
        green_channel = scale_image(green_channel)
    if alpha_channel is None:
        alpha_channel = np.zeros(red_channel.shape, dtype="uint8") + 255
    else:
        alpha_channel = scale_image(alpha_channel)
    image = np.array([red_channel, green_channel, blue_channel, alpha_channel])
    image = image.transpose().copy()  # Have to make sure it's contiguous
    pw.write_png(image, fn)


def write_bitmap(bitmap_array, filename, max_val=None, transpose=False):
    r"""Write out a bitmapped image directly to a PNG file.

    This accepts a three- or four-channel `bitmap_array`.  If the image is not
    already uint8, it will be scaled and converted.  If it is four channel,
    only the first three channels will be scaled, while the fourth channel is
    assumed to be in the range of [0,1]. If it is not four channel, a fourth
    alpha channel will be added and set to fully opaque.  The resultant image
    will be directly written to `filename` as a PNG with no colormap applied.
    `max_val` is a value used if the array is passed in as anything other than
    uint8; it will be the value used for scaling and clipping in the first
    three channels when the array is converted.  Additionally, the minimum is
    assumed to be zero; this makes it primarily suited for the results of
    volume rendered images, rather than misaligned projections.

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
    transpose : boolean, optional
        If transpose is False, we assume that the incoming bitmap_array is such
        that the first element resides in the upper-left corner.  If True, the
        first element will be placed in the lower-left corner.
    """
    if len(bitmap_array.shape) != 3 or bitmap_array.shape[-1] not in (3, 4):
        raise RuntimeError(
            "Expecting image array of shape (N,M,3) or "
            "(N,M,4), received %s" % str(bitmap_array.shape)
        )

    if bitmap_array.dtype != np.uint8:
        s1, s2 = bitmap_array.shape[:2]
        if bitmap_array.shape[-1] == 3:
            alpha_channel = 255 * np.ones((s1, s2, 1), dtype="uint8")
        else:
            alpha_channel = (255 * bitmap_array[:, :, 3]).astype("uint8")
            alpha_channel.shape = s1, s2, 1
        if max_val is None:
            max_val = bitmap_array[:, :, :3].max()
        bitmap_array = np.clip(bitmap_array[:, :, :3] / max_val, 0.0, 1.0) * 255
        bitmap_array = np.concatenate(
            [bitmap_array.astype("uint8"), alpha_channel], axis=-1
        )
    if transpose:
        bitmap_array = bitmap_array.swapaxes(0, 1).copy(order="C")
    if filename is not None:
        pw.write_png(bitmap_array, filename)
    else:
        return pw.write_png_to_string(bitmap_array.copy())
    return bitmap_array


def write_image(image, filename, color_bounds=None, cmap_name=None, func=lambda x: x):
    r"""Write out a floating point array directly to a PNG file, scaling it and
    applying a colormap.

    This function will scale an image and directly call libpng to write out a
    colormapped version of that image.  It is designed for rapid-fire saving of
    image buffers generated using `yt.visualization.api.FixedResolutionBuffers`
    and the likes.

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
        https://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html .
    func : function, optional
        A function to transform the buffer before applying a colormap.

    Returns
    -------
    scaled_image : uint8 image that has been saved

    Examples
    --------

    >>> sl = ds.slice(0, 0.5, "Density")
    >>> frb1 = FixedResolutionBuffer(sl, (0.2, 0.3, 0.4, 0.5), (1024, 1024))
    >>> write_image(frb1[("gas", "density")], "saved.png")
    """
    if cmap_name is None:
        cmap_name = ytcfg.get("yt", "default_colormap")
    if len(image.shape) == 3:
        mylog.info("Using only channel 1 of supplied image")
        image = image[:, :, 0]
    to_plot = apply_colormap(image, color_bounds=color_bounds, cmap_name=cmap_name)
    pw.write_png(to_plot, filename)
    return to_plot


def apply_colormap(image, color_bounds=None, cmap_name=None, func=lambda x: x):
    r"""Apply a colormap to a floating point image, scaling to uint8.

    This function will scale an image and directly call libpng to write out a
    colormapped version of that image.  It is designed for rapid-fire saving of
    image buffers generated using `yt.visualization.api.FixedResolutionBuffers`
    and the likes.

    Parameters
    ----------
    image : array_like
        This is an (unscaled) array of floating point values, shape (N,N,) to
        save in a PNG file.
    color_bounds : tuple of floats, optional
        The min and max to scale between.  Outlying values will be clipped.
    cmap_name : string, optional
        An acceptable colormap.  See either yt.visualization.color_maps or
        https://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html .
    func : function, optional
        A function to transform the buffer before applying a colormap.

    Returns
    -------
    to_plot : uint8 image with colorbar applied.

    """
    if cmap_name is None:
        cmap_name = ytcfg.get("yt", "default_colormap")
    from yt.data_objects.image_array import ImageArray

    image = ImageArray(func(image))
    if color_bounds is None:
        mi = np.nanmin(image[~np.isinf(image)]) * image.uq
        ma = np.nanmax(image[~np.isinf(image)]) * image.uq
        color_bounds = mi, ma
    else:
        color_bounds = [YTQuantity(func(c), image.units) for c in color_bounds]
    image = (image - color_bounds[0]) / (color_bounds[1] - color_bounds[0])
    to_plot = map_to_colors(image, cmap_name)
    to_plot = np.clip(to_plot, 0, 255)
    return to_plot


def map_to_colors(buff, cmap_name):
    lut = get_colormap_lut(cmap_name)

    if isinstance(cmap_name, tuple):
        # If we are using the colorbrewer maps, don't interpolate
        shape = buff.shape
        # We add float_eps so that digitize doesn't go out of bounds
        x = np.mgrid[0.0 : 1.0 + np.finfo(np.float32).eps : lut[0].shape[0] * 1j]
        inds = np.digitize(buff.ravel(), x)
        inds.shape = (shape[0], shape[1])
        mapped = np.dstack([(v[inds] * 255).astype("uint8") for v in lut])
        del inds
    else:
        x = np.mgrid[0.0 : 1.0 : lut[0].shape[0] * 1j]
        mapped = np.dstack([(np.interp(buff, x, v) * 255).astype("uint8") for v in lut])
    return mapped.copy("C")


def strip_colormap_data(
    fn="color_map_data.py",
    cmaps=(
        "jet",
        "cmyt.algae",
        "hot",
        "gist_stern",
        "RdBu",
        "cmyt.pastel",
        "cmyt.kelp",
        "cmyt.arbre",
        "cmyt.octarine",
        "cmyt.dusk",
    ),
):
    import pprint

    from yt._maintenance.deprecation import issue_deprecation_warning

    from . import color_maps as rcm

    issue_deprecation_warning(
        "yt.visualization.image_writer.strip_colormap_data is deprecated.",
        since="4.1.0",
        removal="4.2.0",
    )
    f = open(fn, "w")
    f.write("### Auto-generated colormap tables, taken from Matplotlib ###\n\n")
    f.write("from numpy import array\n")
    f.write("color_map_luts = {}\n\n\n")
    if cmaps is None:
        cmaps = rcm.ColorMaps
    if isinstance(cmaps, str):
        cmaps = [cmaps]
    for cmap_name in sorted(cmaps):
        vals = get_colormap_lut(cmap_name)
        f.write(f"### {cmap_name} ###\n\n")
        f.write(f"color_map_luts['{cmap_name}'] = \\\n")
        f.write("   (\n")
        for v in vals:
            f.write(pprint.pformat(v, indent=3))
            f.write(",\n")
        f.write("   )\n\n")
    f.close()


def splat_points(image, points_x, points_y, contribution=None, transposed=False):
    if contribution is None:
        contribution = 100.0
    val = contribution * 1.0 / points_x.size
    if transposed:
        points_y = 1.0 - points_y
        points_x = 1.0 - points_x
    im = image.copy()
    au.add_points_to_image(im, points_x, points_y, val)
    return im


def write_projection(
    data,
    filename,
    colorbar=True,
    colorbar_label=None,
    title=None,
    vmin=None,
    vmax=None,
    take_log=True,
    figsize=(8, 6),
    dpi=100,
    cmap_name=None,
    extent=None,
    xlabel=None,
    ylabel=None,
):
    r"""Write a projection or volume rendering to disk with a variety of
    pretty parameters such as limits, title, colorbar, etc.  write_projection
    uses the standard matplotlib interface to create the figure.  N.B. This code
    only works *after* you have created the projection using the standard
    framework (i.e. the Camera interface or off_axis_projection).

    Accepts an NxM sized array representing the projection itself as well
    as the filename to which you will save this figure.  Note that the final
    resolution of your image will be a product of dpi/100 * figsize.

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
    vmin : float or None
        the lower limit of the zaxis (part of matplotlib api)
    vmax : float or None
        the lower limit of the zaxis (part of matplotlib api)
    take_log : boolean
        plot the log of the data array (and take the log of the limits if set)?
    figsize : array_like
        width, height in inches of final image
    dpi : int
        final image resolution in pixels / inch
    cmap_name : string
        The name of the colormap.

    Examples
    --------

    >>> image = off_axis_projection(ds, c, L, W, N, "Density", no_ghost=False)
    >>> write_projection(
    ...     image,
    ...     "test.png",
    ...     colorbar_label="Column Density (cm$^{-2}$)",
    ...     title="Offaxis Projection",
    ...     vmin=1e-5,
    ...     vmax=1e-3,
    ...     take_log=True,
    ... )
    """
    if cmap_name is None:
        cmap_name = ytcfg.get("yt", "default_colormap")
    import matplotlib.colors
    import matplotlib.figure

    # If this is rendered as log, then apply now.
    if take_log:
        norm_cls = matplotlib.colors.LogNorm
    else:
        norm_cls = matplotlib.colors.Normalize
    norm = norm_cls(vmin=vmin, vmax=vmax)

    # Create the figure and paint the data on
    fig = matplotlib.figure.Figure(figsize=figsize)
    ax = fig.add_subplot(111)

    cax = ax.imshow(
        data.to_ndarray(),
        norm=norm,
        extent=extent,
        cmap=cmap_name,
    )

    if title:
        ax.set_title(title)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    # Suppress the x and y pixel counts
    if extent is None:
        ax.set_xticks(())
        ax.set_yticks(())

    # Add a color bar and label if requested
    if colorbar:
        cbar = fig.colorbar(cax)
        if colorbar_label:
            cbar.ax.set_ylabel(colorbar_label)

    filename = validate_image_name(filename)
    canvas = get_canvas(fig, filename)

    mylog.info("Saving plot %s", filename)
    fig.tight_layout()

    canvas.print_figure(filename, dpi=dpi)
    return filename


def display_in_notebook(image, max_val=None):
    """
    A helper function to display images in an IPython notebook

    Must be run from within an IPython notebook, or else it will raise
    a YTNotInsideNotebook exception.

    Parameters
    ----------
    image : array_like
        This is an (unscaled) array of floating point values, shape (N,N,3) or
        (N,N,4) to display in the notebook. The first three channels will be
        scaled automatically.
    max_val : float, optional
        The upper limit to clip values of the image.  Only applies to the first
        three channels.
    """

    if "__IPYTHON__" in dir(builtins):
        from IPython.core.displaypub import publish_display_data

        data = write_bitmap(image, None, max_val=max_val)
        publish_display_data(
            data={"image/png": data},
            source="yt.visualization.image_writer.display_in_notebook",
        )
    else:
        raise YTNotInsideNotebook

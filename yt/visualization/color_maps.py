from typing import Union

import cmyt  # noqa: F401
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from yt.funcs import get_brewer_cmap
from yt.utilities.logger import ytLogger as mylog

from . import _colormap_data as _cm

yt_colormaps = {}


def add_colormap(name, cdict):
    """
    Adds a colormap to the colormaps available in yt for this session
    """
    # Note: this function modifies the global variable 'yt_colormaps'
    yt_colormaps[name] = LinearSegmentedColormap(name, cdict, 256)
    mpl.colormaps.register(yt_colormaps[name])


# YTEP-0040 backward compatibility layer
# yt colormaps used to be defined here, but were migrated to an external
# package, cmyt. In the process, 5 of them were renamed. We register them again here
# under their historical names to preserves backwards compatibility.

_HISTORICAL_ALIASES = {
    "arbre": "cmyt.arbre",
    "algae": "cmyt.algae",
    "bds_highcontrast": "cmyt.algae",
    "octarine": "cmyt.octarine",
    "dusk": "cmyt.dusk",
    "kamae": "cmyt.pastel",
    "kelp": "cmyt.kelp",
    "black_blueish": "cmyt.pixel_blue",
    "black_green": "cmyt.pixel_green",
    "purple_mm": "cmyt.xray",
}


def register_yt_colormaps_from_cmyt():
    """
    For backwards compatibility, register yt colormaps without the "cmyt."
    prefix, but do it in a collision-safe way.
    """

    for hist_name, alias in _HISTORICAL_ALIASES.items():
        # note that mpl.colormaps.__getitem__ returns *copies*
        cmap = mpl.colormaps[alias]
        cmap.name = hist_name
        mpl.colormaps.register(cmap)
        mpl.colormaps.register(cmap.reversed())


register_yt_colormaps_from_cmyt()

# Add colormaps in _colormap_data.py that weren't defined here
_vs = np.linspace(0, 1, 256)
for k, v in list(_cm.color_map_luts.items()):
    if k in yt_colormaps:
        continue
    cdict = {
        "red": np.transpose([_vs, v[0], v[0]]),
        "green": np.transpose([_vs, v[1], v[1]]),
        "blue": np.transpose([_vs, v[2], v[2]]),
    }
    try:
        add_colormap(k, cdict)
    except ValueError:
        # expected if another map with identical name was already registered
        mylog.warning("cannot register colormap '%s' (naming collision)", k)


def get_colormap_lut(cmap_id: Union[tuple[str, str], str]):
    # "lut" stands for "lookup table". This function provides a consistent and
    # reusable accessor to a hidden (and by default, uninitialized) attribute
    # (`_lut`) in registered colormaps, from matplotlib or palettable.
    # colormap "lookup tables" are RGBA arrays in matplotlib,
    # and contain sufficient data to reconstruct the colormaps entirely.
    # This exists mostly for historical reasons, hence the custom output format.
    # It isn't meant as part of yt's public api.

    if isinstance(cmap_id, tuple) and len(cmap_id) == 2:
        cmap = get_brewer_cmap(cmap_id)
    elif isinstance(cmap_id, str):
        cmap = mpl.colormaps[cmap_id]
    else:
        raise TypeError(
            "Expected a string or a 2-tuple of strings as a colormap id. "
            f"Received: {cmap_id}"
        )
    if not cmap._isinit:
        cmap._init()
    r = cmap._lut[:-3, 0]
    g = cmap._lut[:-3, 1]
    b = cmap._lut[:-3, 2]
    a = np.ones(b.shape)
    return [r, g, b, a]


def show_colormaps(subset="all", filename=None):
    """
    Displays the colormaps available to yt.  Note, most functions can use
    both the matplotlib and the native yt colormaps; however, there are
    some special functions existing within image_writer.py (e.g. write_image()
    write_bitmap(), etc.), which cannot access the matplotlib
    colormaps.

    In addition to the colormaps listed, one can access the reverse of each
    colormap by appending a "_r" to any map.

    If you wish to only see certain colormaps, include them in the cmap_list
    attribute.

    Parameters
    ----------

    subset : string, or list of strings, optional

        valid values : "all", "yt_native", or list of cmap names
        default : "all"

        As mentioned above, a few functions can only access yt_native
        colormaps.  To display only the yt_native colormaps, set this
        to "yt_native".

        If you wish to only see a few colormaps side by side, you can
        include them as a list of colormap names.
        Example: ['cmyt.algae', 'gist_stern', 'cmyt.kamae', 'nipy_spectral']

    filename : string, opt

        default: None

        If filename is set, then it will save the colormaps to an output
        file.  If it is not set, it will "show" the result interactively.
    """
    from matplotlib import pyplot as plt

    a = np.outer(np.arange(0, 1, 0.01), np.ones(10))
    if subset == "all":
        maps = [
            m
            for m in plt.colormaps()
            if (not m.startswith("idl")) & (not m.endswith("_r"))
        ]
    elif subset == "yt_native":
        maps = [
            m
            for m in _cm.color_map_luts
            if (not m.startswith("idl")) & (not m.endswith("_r"))
        ]
    else:
        try:
            maps = [m for m in plt.colormaps() if m in subset]
            if len(maps) == 0:
                raise AttributeError
        except AttributeError as e:
            raise AttributeError(
                "show_colormaps requires subset attribute "
                "to be 'all', 'yt_native', or a list of "
                "valid colormap names."
            ) from e

    maps = sorted(set(maps))
    # scale the image size by the number of cmaps
    plt.figure(figsize=(2.0 * len(maps) / 10.0, 6))
    plt.subplots_adjust(top=0.7, bottom=0.05, left=0.01, right=0.99)
    l = len(maps) + 1
    for i, m in enumerate(maps):
        plt.subplot(1, l, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect="auto", cmap=mpl.colormaps[m], origin="lower")
        plt.title(m, rotation=90, fontsize=10, verticalalignment="bottom")
    if filename is not None:
        plt.savefig(filename, dpi=100, facecolor="gray")
    else:
        plt.show()


def make_colormap(ctuple_list, name=None, interpolate=True):
    """
    This generates a custom colormap based on the colors and spacings you
    provide.  Enter a ctuple_list, which consists of tuples of (color, spacing)
    to return a colormap appropriate for use in yt.  If you specify a
    name, it will automatically be added to the current session as a valid
    colormap.

    Output colormap is in the format yt expects for adding a colormap to the
    current session: a dictionary with the appropriate RGB channels each
    consisting of a 256x3 array :
    First number is the number at which we are defining a color breakpoint
    Second number is the (0..1) number to interpolate to when coming *from below*
    Third number is the (0..1) number to interpolate to when coming *from above*

    Parameters
    ----------

    ctuple_list: list of (color, float) tuples
        The ctuple_list consists of pairs of (color, interval) tuples
        identifying the colors to use in the colormap and the intervals
        they take to change to the next color in the list.  A color can
        either be a string of the name of a color, or it can be an array
        of 3 floats, each representing the intensity of R, G, and B on
        a scale of 0 to 1.  Valid color names and their equivalent
        arrays are listed below.

        Any interval can be given for the different color tuples, and
        the total of all the intervals will be scaled to the 256 output
        elements.

        If a ctuple_list ends with a color and a non-zero interval,
        a white 0-interval would be added to the end to finish the
        interpolation.  To avoid finishing with white, specify your own
        zero-interval color at the end.

    name: string, optional
        If you wish this colormap to be added as a valid colormap to the
        current session, specify a name here.  Default: None

    interpolate: boolean
        Designates whether or not the colormap will interpolate between
        the colors provided or just give solid colors across the intervals.
        Default: True

    Preset Color Options
    --------------------

    'white' : np.array([255, 255, 255 ])/255.
    'gray' : np.array([130, 130, 130])/255.
    'dgray' : np.array([80, 80, 80])/255.
    'black' : np.array([0, 0, 0])/255.
    'blue' : np.array([0, 0, 255])/255.
    'dblue' : np.array([0, 0, 160])/255.
    'purple' : np.array([100, 0, 200])/255.
    'dpurple' : np.array([66, 0, 133])/255.
    'dred' : np.array([160, 0, 0])/255.
    'red' : np.array([255, 0, 0])/255.
    'orange' : np.array([255, 128, 0])/255.
    'dorange' : np.array([200,100, 0])/255.
    'yellow' : np.array([255, 255, 0])/255.
    'dyellow' : np.array([200, 200, 0])/255.
    'green' : np.array([0, 255, 0])/255.
    'dgreen' : np.array([0, 160, 0])/255.

    Examples
    --------

    To obtain a colormap that starts at black with equal intervals in green,
    blue, red, yellow in that order and interpolation between those colors.
    (In reality, it starts at black, takes an interval of 10 to interpolate to
    green, then an interval of 10 to interpolate to blue, then an interval of
    10 to interpolate to red.)

    >>> cm = make_colormap([("black", 10), ("green", 10), ("blue", 10), ("red", 0)])

    To add a colormap that has five equal blocks of solid major colors to
    the current session as "steps":

    >>> make_colormap(
    ...     [("red", 10), ("orange", 10), ("yellow", 10), ("green", 10), ("blue", 10)],
    ...     name="steps",
    ...     interpolate=False,
    ... )

    To add a colormap that looks like the French flag (i.e. equal bands of
    blue, white, and red) using your own RGB keys, then to display it:

    >>> make_colormap(
    ...     [([0, 0, 1], 10), ([1, 1, 1], 10), ([1, 0, 0], 10)],
    ...     name="french_flag",
    ...     interpolate=False,
    ... )
    >>> show_colormaps(["french_flag"])

    """
    # aliases for different colors
    color_dict = {
        "white": np.array([255, 255, 255]) / 255.0,
        "gray": np.array([130, 130, 130]) / 255.0,
        "dgray": np.array([80, 80, 80]) / 255.0,
        "black": np.array([0, 0, 0]) / 255.0,
        "blue": np.array([0, 0, 255]) / 255.0,
        "dblue": np.array([0, 0, 160]) / 255.0,
        "purple": np.array([100, 0, 200]) / 255.0,
        "dpurple": np.array([66, 0, 133]) / 255.0,
        "dred": np.array([160, 0, 0]) / 255.0,
        "red": np.array([255, 0, 0]) / 255.0,
        "orange": np.array([255, 128, 0]) / 255.0,
        "dorange": np.array([200, 100, 0]) / 255.0,
        "yellow": np.array([255, 255, 0]) / 255.0,
        "dyellow": np.array([200, 200, 0]) / 255.0,
        "green": np.array([0, 255, 0]) / 255.0,
        "dgreen": np.array([0, 160, 0]) / 255.0,
    }

    cmap = np.zeros((256, 3))

    # If the user provides a list with a non-zero final interval, it
    # doesn't make sense because you have an interval but no final
    # color to which it interpolates.  So provide a 0-length white final
    # interval to end the previous interval in white.
    if ctuple_list[-1][1] != 0:
        ctuple_list.append(("white", 0))

    # Figure out how many intervals there are total.
    rolling_index = 0
    for i, (color, interval) in enumerate(ctuple_list):
        if isinstance(color, str):
            ctuple_list[i] = (color_dict[color], interval)
        rolling_index += interval
    scale = 256.0 / rolling_index
    n = len(ctuple_list)

    # Step through each ctuple and interpolate from one color to the
    # next over the interval provided
    rolling_index = 0
    for i in range(n - 1):
        color, interval = ctuple_list[i]
        interval *= scale
        next_index = rolling_index + interval
        next_color, next_interval = ctuple_list[i + 1]

        if not interpolate:
            next_color = color

        # Interpolate the R, G, and B channels from one color to the next
        # Use np.round to make sure you're on a discrete index
        interval = int(np.round(next_index) - np.round(rolling_index))
        for j in np.arange(3):
            cmap[
                int(np.rint(rolling_index)) : int(np.rint(next_index)), j
            ] = np.linspace(color[j], next_color[j], num=interval)

        rolling_index = next_index

    # Return a dictionary with the appropriate RGB channels each consisting of
    # a 256x3 array in the format that is expected by add_colormap() to add a
    # colormap to the session.

    # The format is as follows:
    #   First number is the number at which we are defining a color breakpoint
    #   Second number is the (0..1) number to interpolate to when coming *from below*
    #   Third number is the (0..1) number to interpolate to when coming *from above*
    _vs = np.linspace(0, 1, 256)
    cdict = {
        "red": np.transpose([_vs, cmap[:, 0], cmap[:, 0]]),
        "green": np.transpose([_vs, cmap[:, 1], cmap[:, 1]]),
        "blue": np.transpose([_vs, cmap[:, 2], cmap[:, 2]]),
    }

    if name is not None:
        add_colormap(name, cdict)

    return cdict

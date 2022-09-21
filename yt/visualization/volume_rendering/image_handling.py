import numpy as np

from yt.funcs import mylog
from yt.utilities.on_demand_imports import _h5py as h5py


def export_rgba(
    image,
    fn,
    h5=True,
    fits=False,
):
    """
    This function accepts an *image*, of shape (N,M,4) corresponding to r,g,b,a,
    and saves to *fn*.  If *h5* is True, then it will save in hdf5 format.  If
    *fits* is True, it will save in fits format.
    """
    if (not h5 and not fits) or (h5 and fits):
        raise ValueError("Choose either HDF5 or FITS format!")
    if h5:
        f = h5py.File(f"{fn}.h5", mode="w")
        f.create_dataset("R", data=image[:, :, 0])
        f.create_dataset("G", data=image[:, :, 1])
        f.create_dataset("B", data=image[:, :, 2])
        f.create_dataset("A", data=image[:, :, 3])
        f.close()
    if fits:
        from yt.visualization.fits_image import FITSImageData

        data = {}
        data["r"] = image[:, :, 0]
        data["g"] = image[:, :, 1]
        data["b"] = image[:, :, 2]
        data["a"] = image[:, :, 3]
        fib = FITSImageData(data)
        fib.writeto(f"{fn}.fits", overwrite=True)


def import_rgba(name, h5=True):
    """
    This function will read back in an HDF5 file, as saved by export_rgba, and
    return the frames to the user.  *name* is the name of the file to be read
    in.
    """
    if h5:
        f = h5py.File(name, mode="r")
        r = f["R"].value
        g = f["G"].value
        b = f["B"].value
        a = f["A"].value
        f.close()
    else:
        mylog.error("No support for fits import.")
    return np.array([r, g, b, a]).swapaxes(0, 2).swapaxes(0, 1)


def plot_channel(
    image,
    name,
    cmap="gist_heat",
    log=True,
    dex=3,
    zero_factor=1.0e-10,
    label=None,
    label_color="w",
    label_size="large",
):
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
    from matplotlib import pyplot as plt
    from matplotlib.colors import LogNorm

    from yt.visualization.color_maps import _get_cmap

    Nvec = image.shape[0]
    image[np.isnan(image)] = 0.0
    ma = image[image > 0.0].max()
    image[image == 0.0] = ma * zero_factor
    if log:
        mynorm = LogNorm(ma / (10.0**dex), ma)

    fig = plt.gcf()
    ax = plt.gca()
    fig.clf()
    fig.set_dpi(100)
    fig.set_size_inches((Nvec / 100.0, Nvec / 100.0))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0, hspace=0.0
    )
    mycm = _get_cmap(cmap)
    if log:
        ax.imshow(image, cmap=mycm, norm=mynorm, interpolation="nearest")
    else:
        ax.imshow(image, cmap=mycm, interpolation="nearest")
    if label is not None:
        ax.text(20, 20, label, color=label_color, size=label_size)
    fig.savefig(f"{name}_{cmap}.png")
    fig.clf()


def plot_rgb(image, name, label=None, label_color="w", label_size="large"):
    """
    This will plot the r,g,b channels of an *image* of shape (N,M,3) or
    (N,M,4).  *name* is the prefix of the file name, which will be supplemented
    with "_rgb.png."  *label*, *label_color* and *label_size* may also be
    specified.
    """
    import matplotlib.pyplot as plt

    Nvec = image.shape[0]
    image[np.isnan(image)] = 0.0
    if image.shape[2] >= 4:
        image = image[:, :, :3]

    fig = plt.gcf()
    ax = plt.gca()
    fig.clf()
    fig.set_dpi(100)
    fig.set_size_inches((Nvec / 100.0, Nvec / 100.0))
    fig.subplots_adjust(
        left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0, hspace=0.0
    )
    ax.imshow(image, interpolation="nearest")
    if label is not None:
        ax.text(20, 20, label, color=label_color, size=label_size)
    fig.savefig(f"{name}_rgb.png")
    fig.clf()

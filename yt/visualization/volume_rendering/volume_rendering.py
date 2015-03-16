"""
Volume rendering

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from scene import Scene
from camera import Camera
from render_source import VolumeSource
from .utils import data_source_or_all
from yt.funcs import mylog


def volume_render(data_source, field=None, fname=None, clip_ratio=None):
    r""" Create a simple volume rendering of a data source.

    A helper function that creates a default camera view, transfer
    function, and image size. Using these, it returns an image and
    an instance of the Scene class, allowing one to further modify
    their rendering.

    Parameters
    ----------
    data_source : :class:`yt.data_objects.data_containers.AMR3DData`
        This is the source to be rendered, which can be any arbitrary yt
        3D object
    field: string, tuple, optional
        The field to be rendered. By default, this will use the first
        field in data_source.ds.field_list.  A default transfer function
        will be built that spans the range of values for that given field,
        and the field will be logarithmically scaled if the field_info
        object specifies as such.
    fname: string, optional
        If specified, the resulting rendering will be saved to this filename
        in png format.
    clip_ratio: float, optional
        If specified, the resulting image will be clipped before saving,
        using a threshold based on clip_ratio multiplied by the standard
        deviation of the pixel values. Recommended values are between 2 and 6.

    Returns
    -------
    im: ImageArray
        The resulting image, stored as an ImageArray object.
    sc: Scene
        A :class:`yt.visualization.volume_rendering.scene.Scene` object
        that was constructed during the rendering. Useful for further
        modifications, rotations, etc.

    Example:
    >>> import yt
    >>> ds = yt.load("Enzo_64/DD0046/DD0046")
    >>> im, sc = yt.volume_render(ds, fname='test.png', clip_ratio=4.0)
    """
    data_source = data_source_or_all(data_source)
    sc = Scene()
    if field is None:
        data_source.ds.index
        for ftype, f in sorted(data_source.ds.field_list):
            if ftype == "all": continue
            print ftype, f
            if f == 'Density': field = (ftype, f)
            elif f == 'density': field = (ftype, f)
            elif ftype != 'index' and not 'particle' in f:
                field = (ftype, f)
                break
        else:
            raise RuntimeError("Could not find default field. Please set explicitly in volume_render call")
        mylog.info('Setting default field to %s' % field.__repr__())

    vol = VolumeSource(data_source, field=field)
    cam = Camera(data_source)
    sc.set_camera(cam)
    sc.add_source(vol)
    im = sc.render(fname=fname, clip_ratio=clip_ratio)
    return im, sc

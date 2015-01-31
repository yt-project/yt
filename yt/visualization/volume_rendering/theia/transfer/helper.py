#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


"""Transfer Function Helper Module
This module supports transfer functions that are defined in this directory.	 General
functions should be placed here that might be utilized by multiple transfer functions.
"""

from yt.visualization.volume_rendering.transfer_functions import ColorTransferFunction


def yt_to_rgba_transfer(yt_tf):
	"""Takes a YT Transfer function (ColorTransferFunction) and forms a NumPy Array.
	This is commonly used to hack together a contiguous array from a Yt
	transfer function.
	"""
	return (yt_tf.red.y, yt_tf.green.y, yt_tf.blue.y, yt_tf.alpha.y)

def bounds_to_scale(bounds):
	"""Defines a transfer function scale factor given a transfer function.
	"""
	mi, ma = bounds

	return (1.0)/(ma - mi)

def bounds_to_offset(bounds):
	"""Defines a transfer function offset given a transfer function.
	"""
	mi, ma = bounds

	return mi

def make_yt_transfer(bounds = (0.0, 1.0), colormap = "algae", bins = 1000, scale = 1.0, scale_func = None):
	"""Defines a transfer function offset given a transfer function.
	"""
	mi, ma = bounds
	transfer = ColorTransferFunction( (mi, ma), bins)
        if scale_func == None :
	    transfer.map_to_colormap(mi, ma, colormap=colormap, scale=1.0)
        else :
	    transfer.map_to_colormap(mi, ma, colormap=colormap, scale_func = scale_func)
            

	return transfer

def base_directory(file):

	dir = os.path.dirname(file)
	src = os.path.join(dir, src)
	return src

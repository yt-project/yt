# -----------------------------------------------------------------------------
# Copyright (c) 2016, Fabian Koller (HZDR)
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np

from yt.utilities.logger import ytLogger as mylog


def parse_unit_dimension(unit_dimension):
    """Transforms an openPMD unitDimension into a string.

    Parameters
    ----------
    unit_dimension : array_like
        integer array of length 7 with one entry for the dimensional component of every SI unit

        [0] length L,
        [1] mass M,
        [2] time T,
        [3] electric current I,
        [4] thermodynamic temperature theta,
        [5] amount of substance N,
        [6] luminous intensity J

    References
    ----------
    .. https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md#unit-systems-and-dimensionality

    Returns
    -------
    str

    Examples
    --------
    >>> velocity = [1., 0., -1., 0., 0., 0., 0.]
    >>> print parse_unit_dimension(velocity)
    'm**1*s**-1'

    >>> magnetic_field = [0., 1., -2., -1., 0., 0., 0.]
    >>> print parse_unit_dimension(magnetic_field)
    'kg**1*s**-2*A**-1'
    """
    if len(unit_dimension) is not 7:
        mylog.error("open_pmd - SI must have 7 base dimensions!")
    unit_dimension = np.asarray(unit_dimension, dtype=np.int)
    dim = []
    si = ["m",
          "kg",
          "s",
          "A",
          "C",
          "mol",
          "cd"]
    for i in range(7):
        if unit_dimension[i] != 0:
            dim.append("{}**{}".format(si[i], unit_dimension[i]))
    return "*".join(dim)


def is_const_component(record_component):
    """Determines whether a group or dataset in the HDF5 file is constant.

    Parameters
    ----------
    record_component : h5py.Group or h5py.Dataset

    Returns
    -------
    bool
        True if constant, False otherwise

    References
    ----------
    .. https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md,
       section 'Constant Record Components'
    """
    return "value" in list(record_component.attrs.keys())

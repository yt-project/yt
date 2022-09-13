import numpy as np

from yt.utilities.logger import ytLogger as mylog


def parse_unit_dimension(unit_dimension):
    r"""Transforms an openPMD unitDimension into a string.

    Parameters
    ----------
    unit_dimension : array_like
        integer array of length 7 with one entry for the dimensional component of every
        SI unit

        [0] length L,
        [1] mass M,
        [2] time T,
        [3] electric current I,
        [4] thermodynamic temperature theta,
        [5] amount of substance N,
        [6] luminous intensity J

    References
    ----------

    https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md#unit-systems-and-dimensionality


    Returns
    -------
    str

    Examples
    --------
    >>> velocity = [1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    >>> print(parse_unit_dimension(velocity))
    'm**1*s**-1'

    >>> magnetic_field = [0.0, 1.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    >>> print(parse_unit_dimension(magnetic_field))
    'kg**1*s**-2*A**-1'
    """
    if len(unit_dimension) != 7:
        mylog.error("SI must have 7 base dimensions!")
    unit_dimension = np.asarray(unit_dimension, dtype="int64")
    dim = []
    si = ["m", "kg", "s", "A", "C", "mol", "cd"]
    for i in np.arange(7):
        if unit_dimension[i] != 0:
            dim.append(f"{si[i]}**{unit_dimension[i]}")
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
    return "value" in record_component.attrs.keys()


def get_component(group, component_name, index=0, offset=None):
    """Grabs a dataset component from a group as a whole or sliced.

    Parameters
    ----------
    group : h5py.Group
    component_name : str
        relative path of the component in the group
    index : int, optional
        first entry along the first axis to read
    offset : int, optional
        number of entries to read
        if not supplied, every entry after index is returned

    Notes
    -----
    This scales every entry of the component with the respective "unitSI".

    Returns
    -------
    ndarray
        (N,) 1D in case of particle data
        (O,P,Q) 1D/2D/3D in case of mesh data
    """
    record_component = group[component_name]
    unit_si = record_component.attrs["unitSI"]
    if is_const_component(record_component):
        shape = np.asarray(record_component.attrs["shape"])
        if offset is None:
            shape[0] -= index
        else:
            shape[0] = offset
        # component is constant, craft an array by hand
        return np.full(shape, record_component.attrs["value"] * unit_si)
    else:
        if offset is not None:
            offset += index
        # component is a dataset, return it (possibly masked)
        return np.multiply(record_component[index:offset], unit_si)

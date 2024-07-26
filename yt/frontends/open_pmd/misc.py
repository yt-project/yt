import numpy as np

from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _openpmd_api as openpmd_api


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
    """Determines whether an iteration record is constant

    Parameters
    ----------
    record_component : openpmd_api.openpmd_api_cxx.Record_component

    Returns
    -------
    bool
        True if constant, False otherwise

    References
    ----------
    .. https://github.com/openPMD/openPMD-standard/blob/latest/STANDARD.md,
       section 'Constant Record Components'
    """
    # return record_component.constant #this doesn't work
    return "value" in record_component.attributes


def get_coordinate(mesh, record_axis):
    """This helper function maps coordinate string to openpmd_api  indices
        This has no use yet but I was hoping to be able to make cylindrical data

    Parameters
    ----------
    mesh: openpmd_api.openpmd_api_cxx.Mesh
    record_axis: a string which specifies which desired axis

    Returns
    -------
    int
        specifying index of axis
    """
    if isinstance(mesh, openpmd_api.io.openpmd_api_cxx.Mesh):
        if "cartesian" in str(mesh.geometry):
            cart_map = {"x": 0, "y": 1, "z": 2}
            return cart_map[record_axis]
        elif "cylindrical" in str(mesh.geometry):
            raise AttributeError
        elif "spherical" in str(mesh.geometry):
            raise AttributeError
        elif "thetaMode" in str(mesh.geometry):
            raise AttributeError


def get_component(record, record_axis, index=0, extent=None):
    """Grabs a Record Component from a Record as a whole or sliced.

    Parameters
    ----------
    record : openpmd_api_cxx.Record
    record_axis : str
        the openpmd_api_cxx.Record_Component string key, not necessarily a physical axis
    index : int, optional
        first entry along the first axis to read
    extent : int, optional
        number of entries to read
        note that the previous frontend named this variable offset,
        which we thinks adds some confusion.
        If not supplied, every entry after index is returned.

    Notes
    -----
    This scales every entry of the component with the respective "unitSI".

    Returns
    -------
    ndarray
        (N,) 1D in case of particle data
        (O,P,Q) 1D/2D/3D in case of mesh data
    """
    record_component = record[record_axis]
    unit_si = record_component.get_attribute("unitSI")
    if is_const_component(record_component):
        shape = np.asarray(record_component.get_attribute("shape"))
        if extent is None:
            shape[0] -= index
        else:
            shape = extent
        # component is constant, craft an array by hand
        registered = record_component.get_attribute("value")
        return np.full(shape, registered * unit_si)
    else:
        if extent is not None:
            extent += index
            # len(record_component.shape) gives us on-disk dimensions of component
            if len(record_component.shape) == 3:
                registered = record_component[
                    index[0] : extent[0], index[1] : extent[1], index[2] : extent[2]
                ]
            elif len(record_component.shape) == 2:
                registered = record_component[
                    index[0] : extent[0], index[1] : extent[1]
                ]
            elif len(record_component.shape) == 1:
                # either a scalar or 1d record component
                registered = record_component[index:extent]
        else:
            # when we don't slice we have to .load_chunk()
            registered = record_component.load_chunk()
        record_component.series_flush()
        return np.multiply(registered, unit_si)

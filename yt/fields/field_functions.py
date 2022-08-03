from inspect import signature
from typing import Callable

import numpy as np

from yt.utilities.lib.misc_utilities import obtain_position_vector


def get_radius(data, field_prefix, ftype):
    unit_system = data.ds.unit_system
    center = data.get_field_parameter("center").in_base(unit_system.name)
    DW = (data.ds.domain_right_edge - data.ds.domain_left_edge).in_base(
        unit_system.name
    )
    # This is in cm**2 so it can be the destination for our r later.
    radius2 = data.ds.arr(
        np.zeros(data[ftype, field_prefix + "x"].shape, dtype="float64"), "cm**2"
    )
    r = radius2.copy()
    if any(data.ds.periodicity):
        rdw = radius2.copy()
    for i, ax in enumerate("xyz"):
        # This will coerce the units, so we don't need to worry that we copied
        # it from a cm**2 array.
        np.subtract(
            data[ftype, f"{field_prefix}{ax}"].in_base(unit_system.name),
            center[i],
            r,
        )
        if data.ds.periodicity[i]:
            np.abs(r, r)
            np.subtract(r, DW[i], rdw)
            np.abs(rdw, rdw)
            np.minimum(r, rdw, r)
        np.multiply(r, r, r)
        np.add(radius2, r, radius2)
        if data.ds.dimensionality < i + 1:
            break
    # Now it's cm.
    np.sqrt(radius2, radius2)
    # Alias it, just for clarity.
    radius = radius2
    return radius


def get_periodic_rvec(data):
    coords = obtain_position_vector(data).d
    if sum(data.ds.periodicity) == 0:
        return coords
    le = data.ds.domain_left_edge.in_units("code_length").d
    dw = data.ds.domain_width.in_units("code_length").d
    for i in range(coords.shape[0]):
        if not data.ds.periodicity[i]:
            continue
        coords[i, ...] -= le[i]
        # figure out which measure is less
        mins = np.argmin(
            [
                np.abs(np.mod(coords[i, ...], dw[i])),
                np.abs(np.mod(coords[i, ...], -dw[i])),
            ],
            axis=0,
        )
        temp_coords = np.mod(coords[i, ...], dw[i])

        # Where second measure is better, updating temporary coords
        ii = mins == 1
        temp_coords[ii] = np.mod(coords[i, ...], -dw[i])[ii]

        # Putting the temporary coords into the actual storage
        coords[i, ...] = temp_coords

        coords[i, ...] + le[i]

    return coords


def validate_field_function(function: Callable) -> None:
    """
    Inspect signature, raise a TypeError if invalid, return None otherwise.
    """
    # This is a helper function to user-intended field registration methods
    # (e.g. Dataset.add_field and yt.derived_field)
    # it is not used in FieldInfoContainer.add_field to optimize performance
    # (inspect.signature is quite expensive and we don't want to validate yt's
    # internal code every time a dataset's fields are defined).

    # lookup parameters that do not have default values
    fparams = signature(function).parameters
    nodefaults = tuple(p.name for p in fparams.values() if p.default is p.empty)
    if nodefaults != ("field", "data"):
        raise TypeError(
            f"Received field function {function} with invalid signature. "
            f"Expected exactly 2 positional parameters ('field', 'data'), got {nodefaults!r}"
        )
    if any(
        fparams[name].kind == fparams[name].KEYWORD_ONLY for name in ("field", "data")
    ):
        raise TypeError(
            f"Received field function {function} with invalid signature. "
            "Parameters 'field' and 'data' must accept positional values "
            "(they cannot be keyword-only)"
        )

import warnings
from collections.abc import Callable
from inspect import Parameter, signature

import numpy as np

from yt.utilities.lib.misc_utilities import obtain_position_vector


def get_radius(data, field_prefix, ftype):
    center = data.get_field_parameter("center").to("code_length")
    DW = (data.ds.domain_right_edge - data.ds.domain_left_edge).to("code_length")
    # This is in code_length so it can be the destination for our r later.
    radius2 = data.ds.arr(
        np.zeros(data[ftype, field_prefix + "x"].shape, dtype="float64"), "code_length"
    )

    r = np.empty_like(radius2, subok=False)
    if any(data.ds.periodicity):
        rdw = radius2.v
    for i, ax in enumerate("xyz"):
        pos = data[ftype, f"{field_prefix}{ax}"]
        if str(pos.units) != "code_length":
            pos = pos.to("code_length")
        np.subtract(
            pos.d,
            center[i].d,
            r,
        )
        if data.ds.periodicity[i]:
            np.abs(r, r)
            np.subtract(r, DW.d[i], rdw)
            np.abs(rdw, rdw)
            np.minimum(r, rdw, r)
        np.multiply(r, r, r)
        np.add(radius2.d, r, radius2.d)
        if data.ds.dimensionality < i + 1:
            break
    # Using the views into the array is not changing units and as such keeps
    # from having to do symbolic manipulations
    np.sqrt(radius2.d, radius2.d)
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
    if "data" not in fparams:
        raise TypeError(
            f"Received field function {function} with invalid signature. "
            f"Expected a 'data' parameter, found {tuple(p.name for p in fparams.values())!r}"
        )

    nodefaults = tuple(p.name for p in fparams.values() if p.default is Parameter.empty)
    if not set(nodefaults).issubset({"field", "data"}):
        raise TypeError(
            f"Received field function {function} with invalid signature. "
            "Expected parameters without default values to be 'data' "
            f"(and optionally, 'field'). Found {nodefaults!r}"
        )

    if any(
        nokeyword := tuple(
            name
            for name in nodefaults
            if fparams[name].kind
            not in (Parameter.KEYWORD_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
        )
    ):
        raise TypeError(
            f"Received field function {function} with invalid signature. "
            f"Expected {nokeyword} to be allowed as a keyword argument."
        )

    defaults = tuple(
        name
        for name in ("data", "field")
        if name in fparams and fparams[name].default is not Parameter.empty
    )
    if defaults:
        warnings.warn(
            f"Received field function {function} with default values for parameters {defaults!r}. "
            "These default values will never be used. Drop them to avoid this warning.",
            UserWarning,
            stacklevel=2,
        )

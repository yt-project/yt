from typing import List

import numpy as np
from numpy.typing import ArrayLike


def _validate_cell_widths(
    cell_widths: List[np.ndarray],
    domain_dimensions: ArrayLike,
) -> List[List[np.ndarray],]:
    # check dimensionality
    ndims = len(np.asarray(domain_dimensions))  # cast to array for type checking
    nwids = len(cell_widths)
    if nwids != ndims:
        raise ValueError(
            f"The number of elements in cell_widths ({nwids}) "
            f"must match the number of dimensions ({ndims})."
        )

    # check the dtypes for each dimension, upcast to float64 if needed
    cast_dims = []
    for idim, cell_wid in enumerate(cell_widths):
        if cell_wid.dtype != np.float64:
            cast_dims.append(idim)

    for idim in cast_dims:
        cell_widths[idim] = cell_widths[idim].astype(np.float64)

    # finally, need to return a list of the cell_widths for each grid object.
    # since there is only a single grid, just wrap it in a list.
    cell_widths_out = [
        cell_widths,
    ]
    return cell_widths_out

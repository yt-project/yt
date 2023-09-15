import numpy as np

from yt._typing import DomainDimensions


def _validate_cell_widths(
    cell_widths: list[np.ndarray],
    domain_dimensions: DomainDimensions,
) -> list[list[np.ndarray]]:
    # check dimensionality
    if (nwids := len(cell_widths)) != (ndims := len(domain_dimensions)):
        raise ValueError(
            f"The number of elements in cell_widths ({nwids}) "
            f"must match the number of dimensions ({ndims})."
        )

    # check the dtypes for each dimension, upcast to float64 if needed
    for idim in range(len(cell_widths)):
        cell_widths[idim] = cell_widths[idim].astype(np.float64, copy=False)

    # finally, need to return a list of the cell_widths for each grid object.
    # since there is only a single grid, just wrap it in a list.
    return [cell_widths]

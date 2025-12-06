import numpy as np
import numpy.typing as nptype

from yt._typing import DomainDimensions


def _validate_cell_widths(
    cell_widths: list[nptype.NDArray],
    domain_dimensions: DomainDimensions,
) -> list[nptype.NDArray]:
    # check dimensionality
    if (nwids := len(cell_widths)) != (ndims := len(domain_dimensions)):
        raise ValueError(
            f"The number of elements in cell_widths ({nwids}) "
            f"must match the number of dimensions ({ndims})."
        )

    # check the dtypes for each dimension, upcast to float64 if needed
    for idim in range(len(cell_widths)):
        cell_widths[idim] = cell_widths[idim].astype(np.float64, copy=False)

    return cell_widths

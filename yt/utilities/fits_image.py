from yt.funcs import issue_deprecation_warning
from yt.visualization.fits_image import (
    FITSImageData,
    FITSOffAxisProjection,
    FITSOffAxisSlice,
    FITSProjection,
    FITSSlice,
)

issue_deprecation_warning(
    "The fits_image classes have been moved to yt.visualization "
    "and can be imported from yt directly. yt.utilities.fits_image "
    "is deprecated."
)

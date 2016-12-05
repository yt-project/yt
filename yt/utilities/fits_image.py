from yt.funcs import issue_deprecation_warning

issue_deprecation_warning("The fits_image classes have been moved to yt.visualization "
                          "and are available as top-level imports. yt.utilities.fits_image "
                          "is deprecated.")

from yt.visualization.fits_image import \
    FITSImageData, FITSSlice, FITSProjection, \
    FITSOffAxisSlice, FITSOffAxisProjection
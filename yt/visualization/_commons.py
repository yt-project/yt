import os

from yt.utilities.logger import ytLogger as mylog

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

AGG_FORMATS = [".png", ".jpg", ".jpeg", ".raw", ".rgba", ".tif", ".tiff"]
SUPPORTED_FORMATS = AGG_FORMATS + [".eps", ".ps", ".pdf", ".svg"]


def validate_image_name(filename, suffix: str = ".png") -> str:
    """
    Build a valid image filename with a specified extension (default to png).
    The suffix parameter is ignored if the input filename has a valid extension already.
    Otherwise, suffix is appended to the filename, replacing any existing extension.
    """
    # sanitizing: normalize leading '.'
    suffix = f".{suffix.lstrip('.')}"

    name, psuffix = os.path.splitext(filename)
    if psuffix in SUPPORTED_FORMATS:
        if suffix:
            mylog.warning(
                "Ignoring supplied suffix '%s' in favour of '%s'", suffix, psuffix
            )
        return str(filename)
    elif psuffix:
        mylog.warning(
            "Found unsupported suffix '%s', using '%s' instead.", psuffix, suffix
        )

    return f"{name}{suffix}"


def get_canvas(figure, filename):

    name, suffix = os.path.splitext(filename)

    if not suffix:
        raise ValueError(
            f"Can not determine canvas class from filename '{filename}' "
            f"without an extension."
        )

    if suffix in AGG_FORMATS:
        return FigureCanvasAgg(figure)

    if suffix == ".pdf":
        return FigureCanvasPdf(figure)

    if suffix == ".svg":
        return FigureCanvasSVG(figure)

    if suffix in (".eps", ".ps"):
        return FigureCanvasPS(figure)

    raise ValueError(
        f"No matching canvas for filename '{filename}' with extension '{suffix}'."
    )

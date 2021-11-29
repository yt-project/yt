import os
import warnings
from typing import Optional

import matplotlib
from packaging.version import Version

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

MPL_VERSION = Version(matplotlib.__version__)

DEFAULT_FONT_PROPERTIES = {
    "family": "stixgeneral",
    "size": 18,
}

if MPL_VERSION >= Version("3.4"):
    DEFAULT_FONT_PROPERTIES["math_fontfamily"] = "cm"


AGG_FORMATS = [".png", ".jpg", ".jpeg", ".raw", ".rgba", ".tif", ".tiff"]
SUPPORTED_FORMATS = AGG_FORMATS + [".eps", ".ps", ".pdf", ".svg"]


def normalize_extension_string(s: str) -> str:
    return f".{s.lstrip('.')}"


def validate_image_name(filename, suffix: Optional[str] = None) -> str:
    """
    Build a valid image filename with a specified extension (default to png).
    The suffix parameter is ignored if the input filename has a valid extension already.
    Otherwise, suffix is appended to the filename, replacing any existing extension.
    """
    name, psuffix = os.path.splitext(filename)
    if psuffix in SUPPORTED_FORMATS:
        if suffix is not None:
            suffix = normalize_extension_string(suffix)
        if suffix in SUPPORTED_FORMATS and suffix != psuffix:
            warnings.warn(
                f"Received two valid image formats {psuffix!r} (from filename) "
                f"and {suffix!r} (from suffix). The former is ignored."
            )
            return f"{name}{suffix}"
        return str(filename)

    if suffix is None:
        suffix = ".png"

    suffix = normalize_extension_string(suffix)

    if suffix not in SUPPORTED_FORMATS:
        raise ValueError(f"Unsupported file format {suffix!r}")

    return f"{filename}{suffix}"


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

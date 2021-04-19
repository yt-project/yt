import os
from typing import Optional

from yt.utilities.logger import ytLogger as mylog

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

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
            mylog.warning(
                "Received two valid image formats '%s' (from `filename`) "
                "and '%s' (from `suffix`). The former is ignored.",
                psuffix,
                suffix,
            )
            return f"{name}{suffix}"
        return str(filename)

    if suffix is None:
        suffix = ".png"

    suffix = normalize_extension_string(suffix)

    if suffix not in SUPPORTED_FORMATS:
        raise ValueError("Unsupported file format '{suffix}'.")

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

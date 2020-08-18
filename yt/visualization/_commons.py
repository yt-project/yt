from pathlib import Path

from ._mpl_imports import (
    FigureCanvasAgg,
    FigureCanvasPdf,
    FigureCanvasPS,
    FigureCanvasSVG,
)

AGG_FORMATS = [".png", ".jpg", ".jpeg", ".raw", ".rgba", ".tif", ".tiff"]
SUPPORTED_FORMATS = AGG_FORMATS + [".eps", ".ps", ".pdf", ".svg"]


def validate_image_name(filename, suffix=".png"):
    """
    Build a valid image filename with a specified extension (default to png).
    The suffix parameter is ignored if the input filename has a valid extension already.
    Otherwise, suffix is appended to the filename, replacing any existing extension.
    """
    fn = Path(filename)
    if fn.suffix in SUPPORTED_FORMATS:
        return str(filename)

    if not suffix.startswith("."):
        suffix = f".{suffix}"

    return str(fn.with_suffix(suffix))


def get_canvas(figure, filename):

    suffix = Path(filename).suffix

    if suffix == "":
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
        f"No matching canvas for filename '{filename}' with extension '{suffix}'"
    )

from pathlib import Path

from yt.utilities.logger import ytLogger

from ._mpl_imports import FigureCanvasAgg, FigureCanvasPdf, FigureCanvasPS

SUPPORTED_IMAGE_SUFFIXES = [".png", ".eps", ".ps", ".pdf", ".jpg", ".jpeg"]


def validate_image_name(filename, ext=".png"):
    """
    Build a valid image filename with a specified extension (default to png).
    The ext parameter is ignored if the input filename has a valid extension already.
    Otherwise, ext is appended to filename, taking the place of any existing extension.
    """
    fn = Path(filename)
    if fn.suffix in SUPPORTED_IMAGE_SUFFIXES:
        return str(filename)

    if not ext.startswith("."):
        ext = f".{ext}"

    return str(fn.with_suffix(ext))


def get_canvas(figure, filename, default=None):

    suffix = Path(filename).suffix

    if suffix == ".png":
        return FigureCanvasAgg(figure)

    if suffix == ".pdf":
        return FigureCanvasPdf(figure)

    if suffix in (".eps", ".ps"):
        return FigureCanvasPS(figure)

    if default is not None:
        ytLogger.warning("Unknown suffix %s, using default canvas", suffix)
        return default

    raise ValueError(
        f"No matching canvas for filename {filename} with extension {suffix}"
    )

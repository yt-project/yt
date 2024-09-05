from io import BytesIO

import PIL
from PIL import Image
from PIL.PngImagePlugin import PngInfo

from .._version import __version__ as yt_version


def call_png_write_png(buffer, fileobj, dpi):
    metadata = PngInfo()
    metadata.add_text("Software", f"PIL-{PIL.__version__}|yt-{yt_version}")
    Image.fromarray(buffer).save(
        fileobj, dpi=(dpi, dpi), format="png", pnginfo=metadata
    )


def write_png(buffer, filename, dpi=100):
    with open(filename, "wb") as fileobj:
        call_png_write_png(buffer, fileobj, dpi)


def write_png_to_string(buffer, dpi=100):
    fileobj = BytesIO()
    call_png_write_png(buffer, fileobj, dpi)
    png_str = fileobj.getvalue()
    fileobj.close()
    return png_str

from io import BytesIO

from PIL import Image


def call_png_write_png(buffer, fileobj, dpi):
    Image.fromarray(buffer).save(fileobj, dpi=(dpi, dpi), format="png")


def write_png(buffer, filename, dpi=100):
    with open(filename, "wb") as fileobj:
        call_png_write_png(buffer, fileobj, dpi)


def write_png_to_string(buffer, dpi=100):
    fileobj = BytesIO()
    call_png_write_png(buffer, fileobj, dpi)
    png_str = fileobj.getvalue()
    fileobj.close()
    return png_str

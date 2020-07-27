from io import BytesIO

try:
    # matplotlib switched from an internal submodule _png to using pillow (PIL)
    # between v3.1.0 and v3.3.0
    # So PIL should be available on any system where matplotlib._png doesn't exist
    import matplotlib._png as _png
except ImportError:
    from PIL import Image


def call_png_write_png(buffer, fileobj, dpi):
    try:
        _png.write_png(buffer, fileobj, dpi)
    except NameError:
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

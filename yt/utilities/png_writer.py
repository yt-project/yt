import matplotlib._png as _png

from io import BytesIO

def call_png_write_png(buffer, width, height, fileobj, dpi):
    _png.write_png(buffer, fileobj, dpi)

def write_png(buffer, filename, dpi=100):
    width = buffer.shape[1]
    height = buffer.shape[0]
    with open(filename, "wb") as fileobj:
        call_png_write_png(buffer, width, height, fileobj, dpi)

def write_png_to_string(buffer, dpi=100, gray=0):
    width = buffer.shape[1]
    height = buffer.shape[0]
    fileobj = BytesIO()
    call_png_write_png(buffer, width, height, fileobj, dpi)
    png_str = fileobj.getvalue()
    fileobj.close()
    return png_str

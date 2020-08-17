from pathlib import Path

SUPPORTED_IMAGE_SUFFIXES = [".png", ".eps", ".ps", ".pdf", ".jpg", ".jpeg"]


def validate_image_name(filename, ext=".png"):

    fn = Path(filename)
    if fn.suffix in SUPPORTED_IMAGE_SUFFIXES:
        return str(filename)

    if not ext.startswith("."):
        ext = f".{ext}"

    return str(fn.with_suffix(ext))

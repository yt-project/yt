from io import BytesIO
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image

from yt.utilities.png_writer import write_png, write_png_to_string


def test_write_png_to_string_round_trips_image_and_metadata():
    buffer = np.array(
        [
            [[255, 0, 0], [0, 255, 0]],
            [[0, 0, 255], [255, 255, 255]],
        ],
        dtype="uint8",
    )

    png_bytes = write_png_to_string(buffer, dpi=72)

    assert png_bytes.startswith(b"\x89PNG\r\n\x1a\n")

    with Image.open(BytesIO(png_bytes)) as image:
        assert image.size == (2, 2)
        assert image.info["dpi"] == (72.009, 72.009)
        assert "PIL-" in image.info["Software"]
        assert "|yt-" in image.info["Software"]


def test_write_png_writes_file():
    buffer = np.array([[[1, 2, 3]]], dtype="uint8")
    with TemporaryDirectory() as tmpdir:
        filename = f"{tmpdir}/pixel.png"

        write_png(buffer, filename, dpi=100)

        with Image.open(filename) as image:
            assert image.size == (1, 1)
            assert image.getpixel((0, 0)) == (1, 2, 3)

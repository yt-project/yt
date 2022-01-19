from typing import List, Tuple

import numpy as np


class ZBuffer:
    """A container object for z-buffer arrays

    A zbuffer is a companion array for an image that allows the volume rendering
    infrastructure to determine whether one opaque source is in front of another
    opaque source.  The z buffer encodes the distance to the opaque source
    relative to the camera position.

    Parameters
    ----------
    rgba: MxNx4 image
        The image the z buffer corresponds to
    z: MxN image
        The z depth of each pixel in the image. The shape of the image must be
        the same as each RGBA channel in the original image.

    Examples
    --------
    >>> import numpy as np
    >>> shape = (64, 64)
    >>> b1 = Zbuffer(np.random.random(shape), np.ones(shape))
    >>> b2 = Zbuffer(np.random.random(shape), np.zeros(shape))
    >>> c = b1 + b2
    >>> np.all(c.rgba == b2.rgba)
    True
    >>> np.all(c.z == b2.z)
    True
    >>> np.all(c == b2)
    True

    """

    def __init__(self, rgba, z):
        super().__init__()
        assert rgba.shape[: len(z.shape)] == z.shape
        self.rgba = rgba
        self.z = z
        self.shape = z.shape

    def __add__(self, other):
        assert self.shape == other.shape
        f = self.z < other.z
        if self.z.shape[1] == 1:
            # Non-rectangular
            rgba = self.rgba * f[:, None, :]
            rgba += other.rgba * (1.0 - f)[:, None, :]
        else:
            b = self.z > other.z
            rgba = np.zeros(self.rgba.shape)
            rgba[f] = self.rgba[f]
            rgba[b] = other.rgba[b]
        z = np.min([self.z, other.z], axis=0)
        return ZBuffer(rgba, z)

    def __iadd__(self, other):
        tmp = self + other
        self.rgba = tmp.rgba
        self.z = tmp.z
        return self

    def __eq__(self, other):
        equal = True
        equal *= np.all(self.rgba == other.rgba)
        equal *= np.all(self.z == other.z)
        return equal

    def paint(self, ind, value, z):
        if z < self.z[ind]:
            self.rgba[ind] = value
            self.z[ind] = z


if __name__ == "__main__":
    shape: Tuple[int, ...] = (64, 64)
    shapes: List[Tuple[int, ...]] = [(64, 64), (16, 16, 4), (128,), (16, 32)]
    for shape in shapes:
        b1 = ZBuffer(np.random.random(shape), np.ones(shape))
        b2 = ZBuffer(np.random.random(shape), np.zeros(shape))
        c = b1 + b2
        assert np.all(c.rgba == b2.rgba)
        assert np.all(c.z == b2.z)
        assert np.all(c == b2)

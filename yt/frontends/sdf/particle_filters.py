import numpy as np


def bbox_filter(left, right):

    def myfilter(chunk, mask=None):
        pos = np.array([chunk['x'], chunk['y'], chunk['z']]).T

        # Now get all particles that are within the bbox
        if mask is None:
            mask = np.all(pos >= left, axis=1) * np.all(pos < right, axis=1)
        else:
            np.multiply(mask, np.all(pos >= left, axis=1), mask)
            np.multiply(mask, np.all(pos < right, axis=1), mask)
        return mask

    return myfilter

def sphere_filter(center, radius):

    def myfilter(chunk, mask=None):
        pos = np.array([chunk['x'], chunk['y'], chunk['z']]).T

        # Now get all particles that are within the radius
        if mask is None:
            mask = ((pos-center)**2).sum(axis=1)**0.5 < radius
        else:
            np.multiply(mask, np.linalg.norm(pos - center, 2) < radius, mask)
        return mask

    return myfilter

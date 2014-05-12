import numpy as np


def bbox_filter(left, right, domain_width):

    def myfilter(chunk, mask=None):
        pos = np.array([chunk['x'], chunk['y'], chunk['z']]).T

        # This hurts, but is useful for periodicity. Probably should check first
        # if it is even needed for a given left/right
        for i in range(3):
            pos[:,i] = np.mod(pos[:,i] - left[i], domain_width[i]) + left[i]

        # Now get all particles that are within the bbox
        if mask is None:
            mask = np.all(pos >= left, axis=1) * np.all(pos < right, axis=1)
        else:
            np.multiply(mask, np.all(pos >= left, axis=1), mask)
            np.multiply(mask, np.all(pos < right, axis=1), mask)
        return mask

    return myfilter

def sphere_filter(center, radius, domain_width):

    def myfilter(chunk, mask=None):
        pos = np.array([chunk['x'], chunk['y'], chunk['z']]).T

        # This hurts, but is useful for periodicity. Probably should check first
        # if it is even needed for a given left/right
        for i in range(3):
            pos[:,i] = np.mod(pos[:,i] - left[i], domain_width[i]) + left[i]

        # Now get all particles that are within the radius
        if mask is None:
            mask = ((pos-center)**2).sum(axis=1)**0.5 < radius
        else:
            np.multiply(mask, np.linalg.norm(pos - center, 2) < radius, mask)
        return mask

    return myfilter

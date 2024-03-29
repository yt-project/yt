import numpy as np

# Now define convenience functions


def blankRecordArray(desc, elements):
    """
    Accept a descriptor describing a recordarray, and return one that's full of
    zeros

    This seems like it should be in the numpy distribution...
    """
    blanks = []
    for atype in desc["formats"]:
        blanks.append(np.zeros(elements, dtype=atype))
    return np.rec.fromarrays(blanks, **desc)

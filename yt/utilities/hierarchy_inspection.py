# -*- coding: utf-8 -*-
import inspect

def find_lowest_subclass(candidates):
    """
    This function searches an inheritance hierarchy and looks for the lowest 
    subclass in that hierarchy. If the tree diverges then an exception is raised.
    
    Parameters
    ----------
    candidates : iterable
        An interable object that is a collection of classes to find the lowest
        subclass of.
    
    Returns
    -------
    result : object
        The object in candidates which is the lowest in the tree.
    """
    candidates.sort()
    mros = [set(inspect.getmro(c)) for c in candidates]
    def swipe(mros):
        uniques = []
        for m in mros[1:]:
            x = mros[0].symmetric_difference(m)
            uniques.append(x)
    
        while len(uniques) >1:
            uniques = swipe(uniques)
        return uniques
    
    result = swipe(mros)[0]
    
    if len(result) > 1:
        raise TypeError("A Diverging hierarchy was detected.")
    
    else:
        return list(result)[0]

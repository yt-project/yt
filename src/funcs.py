"""
Useful functions.  If non-original, see function for citation.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

import time
from progressbar import *

def iterable(obj):
    """
    Grabbed from Python Cookbook / matploblib.cbook
    """
    try: len(obj)
    except: return False
    return True

def time_execution(func):
    def wrapper(*arg, **kw):
        t1 = time.time()
        res = func(*arg, **kw)
        t2 = time.time()
        mylog.debug('%s took %0.3f s', func.func_name, (t2-t1))
        return res
    from yt import ytcfg
    from yt import lagosLogger as mylog
    if ytcfg.getboolean("yt","timefunctions") == True:
        return wrapper
    else:
        return func

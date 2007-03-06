"""
Not yet functional.  Better strategy needed first...

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Implement!
@todo: Figure out SQLite/MySQL status
@todo: Unique run ids?
"""

from yt.lagos import *

def submitRun(run):
    """
    This accepts an EnzoRun instance and puts it into the database

    @param run: the EnzoRun instance
    @type run: L{EnzoRun<yt.lagos.EnzoRun>}
    """
    return

def fetchRun(md):
    """
    This takes a metaData string and then returns the EnzoRun instance
    associated with that metadata string

    @param md: metaData
    @type md: string
    """
    return

def submitParameterFile(pf):
    """
    This takes a parameterfile, which should have an associated EnzoRun
    instance, and then it inserts it into the database
    @param pf: parameterfile
    @type pf: L{EnzoHierarchy<EnzoHierarchy>}
    """
    return

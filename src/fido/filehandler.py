"""
High-level file handling takes place here.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Moving to storage
@todo: Update DB with locations of files, as needed
"""

import os, os.path, time, sys, shutil
from glob import glob
from yt.fido import *

def moveFiles(path, basename, extraFiles = [], wd=None):
    orig_wd = os.getcwd()
    print wd
    if wd:
        os.chdir(wd)
    print("ls -1f |grep ^'%s'|grep -v \.dir|xargs mv --target-directory='%s'" \
                       % (basename, path))
    retVal = os.system("ls -1f |grep ^'%s'|grep -v \.dir|xargs mv --target-directory='%s'" \
                       % (basename, path))
    # Log here
    mylog.info("Moved %s to %s", basename, path)
    for file in extraFiles:
        try:
            shutil.copy2(file, path+os.path.sep+file)
        except:
            pass
    os.chdir(orig_wd)
    return retVal

def deleteFiles(path, basename):
    """
    This is the much simpler deletion method, which doesn't do any checks or
    involve any run instantiation, but instead just deletes everything matching the
    pattern.

    @param path: the path where the files currently reside
    @param basename: the basename of the files to delete
    """
    orig_wd = os.getcwd()
    os.chdir(path)
    mylog.info("Deleting everything in %s matching '^%s' (%s)", path, basename,os.getcwd())
    retVal = os.system("ls -1f |grep ^'%s'|grep -v \.dir|xargs rm" % (basename))
    os.chdir(orig_wd)
    mylog.info("Finished with retVal %s", retVal)
    return retVal

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
    """
    This copies select files, based on the basename.  This is used to 'bury' an
    output, and should only be executed when the output exists with other files
    mixed in -- like, for instance, during an enzo run.

    @param path: where to copy the files
    @type path: string
    @param basename: the base name of the output (DataDump0001, etc)
    @type basename: string
    @param extraFiles: extra files to copy (I{not} move)
    @type extraFiles: list of strings
    @param wd: the working directory we should be in to do this
    @type wd: string
    """
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

def moveOutput(source, dest):
    """
    This moves an entire directory tree.  Note that this should be used when an
    output is already 'buried' in a tree.  This will retain the 'buried'
    status.  Should be used as follows:

    moveOutput("/data/my_run1/DataDump0001.dir","/data/my_run2/")

    @param source: the directory tree we're moving
    @type source: string
    @param dest: the place we're moving to
    @type dest: string
    """
    shutil.move(source, os.path.join(os.path.basename(dest,source)))
    mylog.info("Moved %s to %s", basename, path)
    
    

def copyOutput(source, dest):
    """
    This copies a buried output to another buried output location.  Use it
    something like:

    copyOutput("/data/my_run1/DataDump0001.dir","/data/my_run2/")

    @param source: the source location (including the directory name, but not the file names)
    @type source: string
    @param dest: the destination, not including the directory name
    @type dest: string
    """
    bn = os.path.basename(source)
    bnf = bn[:bn.rfind(".")]
    new_location = os.path.join(dest, bn)
    if not os.path.exists(new_location):
        os.makedirs(new_location)
    wd = os.getcwd()
    os.chdir(source)
    #print "ls -1f|grep ^'%s'|xargs cp --target-directory='%s'" % (bnf,new_location)
    os.system("ls -1f|grep ^'%s'|xargs cp --target-directory='%s'" % (bnf,new_location))
    #shutil.copy2(source, new_location)
    os.chdir(wd)
    mylog.info("Copied %s to %s", source, dest)
    return new_location

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

def deleteOutput(outputLocation):
    """
    Very brutal directory deletion command.  Use with care, dudes.

    @param outputLocation: the directory tree to trash
    @type outputLocation: string
    """
    mylog.warning("Okay, as requested, I'm going to delete %s", outputLocation)
    shutil.rmtree(outputLocation)
    mylog.info("Did it!  I hope you're happy.")

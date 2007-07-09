"""
We operate on some very simple principles here.
"""

from yt.fido import *

def moveOutput(basename, newLocation, extraFiles=[]):
    moveGlob("%s*" % (basename), newLocation)
    for file in extraFiles:
        moveGlob(file, location)

def deleteOutput(basename, extraFiles=[]):
    deleteGlob("%s*" % os.path.normpath(filename), location)
    for file in extraFiles:
        deleteGlob(file)

def buryOutput(filename, newLocation=None, extraFiles=[]):
    if not newLocation:
        dirname = "%s.dir" % filename
        newLocation = os.path.abspath(os.path.join(os.getcwd(), dirname))
    moveGlob("%s*" % filename, newLocation)
    for file in extraFiles:
        moveGlob(file, newLocation)
    return os.path.join(newLocation, filename)

def digupOutput(filename, newLocation=None, extraFiles=[]):
    if not newLocation:
        newLocation = getParentDir(filename)
        print newLocation
    else:
        # Make sure we have no relative references
        newLocation = os.path.normpath(newLocation)
    moveGlob("%s*" % filename, newLocation)
    
def moveGlob(globPattern, newLocation):
    if not os.path.isdir(newLocation):
        os.makedirs(newLocation)
    # This is slower than dumping to shell.
    # But, I trust it.  And it should work on boh
    # BSD and GNU util systems.
    for file in glob.glob(globPattern):
        mylog.debug("Moving %s to %s", file, newLocation)
        shutil.move(file, newLocation)

def deleteGlob(globPattern):
    for file in glob.glob(globPattern):
        print "Removing %s" % (file)
        shutil.unlink(file)

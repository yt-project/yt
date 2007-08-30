from yt.fido import *

def copyOutput(basename, newLocation, extraFiles=None):
    if extraFiles == None: extraFiles = []
    copyGlob("%s*" % (basename), newLocation)
    for file in extraFiles:
        copyGlob(file, location)
    return os.path.join(newLocation, os.path.basename(basename))

def moveOutput(basename, newLocation, extraFiles=None):
    if extraFiles == None: extraFiles = []
    moveGlob("%s*" % (basename), newLocation)
    for file in extraFiles:
        moveGlob(file, location)

def deleteOutput(basename, extraFiles=None):
    if extraFiles == None: extraFiles = []
    deleteGlob("%s*" % os.path.normpath(filename), location)
    for file in extraFiles:
        deleteGlob(file)

def buryOutput(filename, newLocation=None, extraFiles=None):
    if extraFiles == None: extraFiles = []
    if not newLocation:
        dirname = "%s.dir" % filename
        newLocation = os.path.abspath(os.path.join(os.getcwd(), dirname))
    moveGlob("%s*" % filename, newLocation)
    for file in extraFiles:
        moveGlob(file, newLocation)
    return os.path.join(newLocation, filename)

def digupOutput(filename, newLocation=None, extraFiles=None):
    if extraFiles == None: extraFiles = []
    if not newLocation:
        newLocation = getParentDir(filename)
    else:
        # Make sure we have no relative references
        newLocation = os.path.normpath(newLocation)
    print "%s*" % (filename), newLocation
    moveGlob("%s*" % filename, newLocation)

def copyGlob(globPattern, newLocation):
    if not os.path.isdir(newLocation):
        os.makedirs(newLocation)
    for file in glob.glob(globPattern):
        mylog.debug("Copying %s to %s", file, newLocation)
        shutil.copy(file, newLocation)
    
def moveGlob(globPattern, newLocation):
    if not os.path.isdir(newLocation):
        os.makedirs(newLocation)
    # This is slower than dumping to shell.
    # But, I trust it.  And it should work on boh
    # BSD and GNU util systems.
    for file in glob.glob(globPattern):
        if os.path.abspath(file) == os.path.abspath(newLocation): continue
        nl = os.path.join(newLocation, os.path.basename(file))
        mylog.debug("Moving %s to %s", file, nl)
        shutil.move(file, nl)

def deleteGlob(globPattern):
    for file in glob.glob(globPattern):
        print "Removing %s" % (file)
        shutil.unlink(file)

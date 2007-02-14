#
#
#

WAITBETWEEN=5
OTHERFILES=["rates.out","cool_rates.out"]

import os, os.path, time, sys, shutil
from glob import glob

def checkForOutput(skipFiles = []):
    newFiles = []
    i = 0
    l = glob("*.hierarchy")
    if len(l) > 0:
        for file in l:
            if file not in skipFiles and os.stat(file).st_size > 0:
                newFiles.append(file.rsplit(".",1)[0])
        # Log
        return newFiles
    return None

def makeNewDirectory(path, basename):
    """
    Makes a new directory, based on the basename and the path

    Arguments
    """
    # Log
    return nn

def moveFiles(path, basename, extraFiles = []):
    retVal = os.system("ls -1f |grep ^'%s'|grep -v \.dir|xargs mv --target-directory='%s'" \
                       % (basename, path))
    # Log here
    for file in extraFiles:
        try:
            shutil.copy2(file, path+os.path.sep+file)
        except:
            pass
    return retVal

def checkForStop():
    if os.path.exists("stopFido"):
        # We should log this rather than print it
        print "FOUND stopFido!"
        os.unlink("stopFido")
        return 1
    return 0

def watchDir(path=".", skip = [], funcHandler = None):
    os.chdir(path)
    doStop = 0
    while not doStop:
        nn = checkForOutput(skip)
        if nn:
            print nn
            for bn in nn:
                print bn
                newDir = "%s.dir" % (bn)
                if not os.path.exists(newDir):
                    os.mkdir(newDir)
                print newDir
                moveFiles(newDir, bn, extraFiles=OTHERFILES)
                if funcHandler:
                    pid=os.fork()
                    if pid == 0:
                        funcHandler(newDir, bn)
                        sys.exit()
        time.sleep(WAITBETWEEN)
        doStop = checkForStop()

if __name__=="__main__":
    watchDir()

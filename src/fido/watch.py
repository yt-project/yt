WAITBETWEEN=5
OTHERFILES=["rates.out","cool_rates.out"]
NEW_OUTPUT_CREATED = "newOutput"

import os, os.path, time, sys, shutil
from glob import glob
from yt.fido import *

def newCheckForOutput(skipFiles = []):
    if os.path.isfile(NEW_OUTPUT_CREATED):
        os.unlink(NEW_OUTPUT_CREATED)
        newFiles = []
        # So something is created!  Now let's snag it
        l = glob("*.hierarchy")
        if len(l) > 0:
            for file in l:
                if file not in skipFiles:
                    newFiles.append(file.rsplit(".",1)[0])
                    mylog.info("Found output %s", newFiles[-1])
            # Do we want to do this?
            return newFiles
    return None

def manualCheckForOutput(skipFiles = []):
    newFiles = []
    i = 0
    l = glob("*.hierarchy")
    if len(l) > 0:
        for file in l:
            if file not in skipFiles and os.stat(file).st_size > 0 \
               and abs(os.stat(file).st_mtime - time.time()) > 5:
                newFiles.append(file.rsplit(".",1)[0])
                mylog.info("Found output %s", newFiles[-1])
        # Log
        return newFiles
    return None

def checkForStop(process=None):
    if os.path.exists("stopFido"):
        # We should log this rather than print it
        mylog.info("Stopping fido")
        os.unlink("stopFido")
        return 1
    if process:
        pp = process.poll()
        if pp != None:
            mylog.info("Process has died; stopping fido")
            return 1
    return 0

def watchDir(path=".", skip = [], funcHandler = None, process=None):
    if process:
        mylog.info("Handed process")
    os.chdir(path)
    doStop = 0
    while not doStop:
        doStop = checkForStop(process)
        nn = newCheckForOutput(skip)
        if nn:
            for bn in nn:
                newDir = "%s.dir" % (bn)
                if not os.path.exists(newDir):
                    os.mkdir(newDir)
                moveFiles(newDir, bn, extraFiles=OTHERFILES)
                if funcHandler:
                    mylog.info("Calling %s", funcHandler.func_name)
                    funcHandler(newDir, bn)
                    mylog.info("Done calling %s, now sleeping", funcHandler.func_name)
        time.sleep(WAITBETWEEN)

if __name__=="__main__":
    watchDir()

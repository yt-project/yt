"""
Output-watcher

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

WAITBETWEEN=5
OTHERFILES=["rates.out","cool_rates.out"]
NEW_OUTPUT_CREATED = "newOutput"

import os, os.path, time, sys, shutil
import yt.lagos as lagos
from yt.fido import *
import glob

def newCheckForOutput(skipFiles = []):
    if os.path.isfile(NEW_OUTPUT_CREATED):
        os.unlink(NEW_OUTPUT_CREATED)
        newFiles = []
        # So something is created!  Now let's snag it
        l = glob.glob("*.hierarchy")
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
    l = glob.glob("*.hierarchy")
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

def watchDir(md=None, path=".", skip = [], funcHandler = None, process=None):
    """
    This is the heart of Fido.  This function watches a given directory for
    output, and when that output is found, it calls moves the output and calls a
    function to deal with the data.

    @keyword path: the path to watch for output
    @keyword skip: skip these files when looking for output
    @type skip: list of strings
    @keyword funcHandler: this function will be called with the EnzoHierarchy instance
    @type funcHandler: function
    @keyword process: watch this process for death (don't use directly)
    @type process: process
    @note: B{This assumes outputs are added in chronological time order!}  This
           function may behave incorrectly if you don't revert and digup runs.
    """
    if process:
        mylog.info("Handed process")
    os.chdir(path)
    if not md:
        md = os.path.basename(os.getcwd())
    thisRun = fetchRun(md, classType=lagos.EnzoParameterFile)
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
                mylog.info("Adding to EnzoRun instance")
                thisRun.addOutputByFilename(os.path.join(newDir,bn))
                submitRun(thisRun)
                if funcHandler:
                    pid = os.fork()
                    if pid:
                        newpid, exit = os.waitpid(pid)
                        mylog.info("Exit status %s from PID %s", exit, newpid)
                    else:
                        mylog.info("Forked process reporting for duty")
                        thisRun.promoteType(-1)
                        mylog.info("Calling %s" % funcHandler.func_name)
                        funcHandler(thisRun.outputs[-1])
                        #thisRun.demoteType(-1)
                        mylog.info("Done calling %s, now dying", funcHandler.func_name)
                        sys.exit()
        time.sleep(WAITBETWEEN)

def runOnce(md=None):
    if not md:
        md = os.path.basename(os.getcwd())
    thisRun = fetchRun(md, classType=lagos.EnzoParameterFile)
    nn = newCheckForOutput()
    for bn in nn:
        newDir = "%s.dir" % (bn)
        if not os.path.exists(newDir):
            os.mkdir(newDir)
        moveFiles(newDir, bn, extraFiles=OTHERFILES)
        mylog.info("Adding %s to EnzoRun instance (%s)", bn, md)
        thisRun.addOutputByFilename(os.path.join(newDir,bn))
        submitRun(thisRun)
    mylog.info("All done, exiting")

if __name__=="__main__":
    watchDir()

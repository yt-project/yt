#
# jobdispatch.py
#
# We define a couple classes and functions to dispatch jobs to hosts.
#

import os, sys, signal

class HostProfile:
    def __init__(self, hostname, submitFunction):
        self.hostname = hostname
        self.submitFunction = submitFunction
    def spawnProcess(self, d):
        self.submitFunction(**d)

def submitLSF(**args):
    """
    Submits to LSF, based on fully optional arguments

    Keyword Arguments:
        queue -- the queue to submit to
        resource -- the resource to submit to
        exe -- the executable to use
    """
    return

def submitDPLACE(wd='.', parameterFile=None, exe="./enzo", restart = False, nproc=1, logFile=None):
    """
    Submits directly to a dplaced mpirun command
    
    Keyword Arguments:
        parameterFile -- the parameter file to initialize with
        exe -- the enzo executable
        restart -- do we need to feed the -r argument?
        nproc -- the number of processors to run on
        logFile -- the logfile
    """
    os.chdir(wd)
    commandLine = "/usr/bin/dplace -s1 /usr/bin/mpirun -np %i" % (nproc)
    commandLine += " %s -d" % (exe)
    if restart:
        commandLine += " -r"
    commandLine += " %s" % (parameterFile)
    if logFile:
        commandLine += " 2>&1 | tee %s" % (logFile)
    print "Executing:"
    print commandLine
    # Additionally, here we will fork out to the watcher process
    # If this is a restart dump, we'll feed a temporary skipFile to the watcher
    # of the restart basename; then after one iteration it'll be taken out of
    # the skipFiles list, and will get moved out
    status = os.system(commandLine)

hostRed = HostProfile("red.slac.stanford.edu", submitDPLACE)

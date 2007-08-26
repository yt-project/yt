"""
Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

This file is part of yt.

yt is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


"""
Profiles for the various hosts we can run on

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Finish testing
"""

from yt.enki import *

import os, sys, signal, subprocess

class HostProfile:
    def __init__(self, hostname, submitFunction):
        self.hostname = hostname
        self.submitFunction = submitFunction
    def spawnProcess(self, d):
        pid = self.submitFunction(**d)
        return pid

def submitLSF(**args):
    """
    Submits to LSF, based on fully optional arguments

    NOT FUNCTIONAL YET.  Will probably use a swig-ified LSF.

    @keyword queue: the queue to submit to
    @keyword resource: the resource to submit to
    @keyword exe: the executable to use
    """
    return

def submitDPLACE(wd='.', parameterFile=None, exe="./enzo", restart = False, nproc=1, logFile=None):
    """
    Submits directly to a dplaced mpirun command
    
    @keyword parameterFile: the parameter file to initialize with
    @keyword exe: the enzo executable
    @keyword restart: do we need to feed the -r argument?
    @keyword nproc: the number of processors to run on
    @type nproc: int
    @keyword logFile: the logfile
    """
    commandLine = "/usr/bin/dplace -s1 /usr/bin/mpirun -np %i" % (nproc)
    commandLine += " %s -d" % (exe)
    if restart:
        commandLine += " -r"
    commandLine += " %s" % (parameterFile)
    if logFile:
        commandLine += " 2>&1 | tee %s" % (logFile)
    mylog.info("Executing: '%s'", commandLine)
    # Additionally, here we will fork out to the watcher process
    # If this is a restart dump, we'll feed a temporary skipFile to the watcher
    # of the restart basename; then after one iteration it'll be taken out of
    # the skipFiles list, and will get moved out
    p = subprocess.Popen(commandLine, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=wd, \
            shell=True, executable="/bin/bash")
    return p

hosts = {}

hosts["red_dplace"] = HostProfile("red.slac.stanford.edu", submitDPLACE)

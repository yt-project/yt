"""
Output-watcher

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
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

from yt.fido import *

class Watcher:
    def __init__(self, title=None, path=".", newPrefix="", oc=None,
                 process=None, functionHandler = None):
        self.originalPath = os.getcwd()
        os.chdir(path)
        if title == None: title = os.path.basename(os.getcwd())
        if oc == None: oc = OutputCollection(title)
        self.title = title
        self.oc = oc
        self.process = process
        self.newPrefix = newPrefix
        self.skipFiles = [] # Forward compatible
        self.functionHandler = None
        if functionHandler != None:
            self.functionHandler = functionHandler()

    def run(self, runOnce=False):
        mylog.info("Entering main Fido loop (CTRL-C or touch 'stopFido' to end)")
        wb = ytcfg.getfloat("fido","WaitBetween")
        while not self.checkForStop():
            nn = self.checkForOutput()
            for bn in nn:
                newName = buryOutput(bn, newPrefix=self.newPrefix)
                self.dealWithOutput(newName)
            if runOnce: break
            try:
                time.sleep(wb)
            except KeyboardInterrupt:
                sys.exit()

    def dealWithOutput(self, filename):
        # First, add it to the OutputCollection
        self.oc.addOutput(filename)
        self.oc.writeOut()
        if self.functionHandler == None:
            return
        # Now, we pass it to our function handler
        pid = os.fork()
        if pid:
            mylog.debug("Waiting on pid %s", pid)
            newpid, exit = os.waitpid(pid,0)
            mylog.debug("Exit status %s from PID %s", exit, newpid)
        else:
            mylog.info("Forked process reporting for duty!")
            self.functionHandler(filename)
            sys.exit()

    def checkForOutput(self):
        newFiles = []
        if os.path.isfile(NEW_OUTPUT_CREATED):
            os.unlink(NEW_OUTPUT_CREATED)
            # So something is created!  Now let's snag it
            # We insert our additional glob patterns here
            filesFound = glob.glob("*.hierarchy") #\
                         #[glob.glob(pat) for pat in GlobPatterns]
            for file in filter(lambda a: a not in self.skipFiles, filesFound):
                #if self.oc.has_key(os.path.basename(file[:-10])): continue
                newFiles.append(file.rsplit(".",1)[0])
                mylog.info("Found output %s", newFiles[-1])
        return newFiles

    def checkForStop(self):
        if os.path.exists("stopFido"):
            # We should log this rather than print it
            mylog.info("Stopping fido")
            os.unlink("stopFido")
            return 1
        if self.process:
            pp = self.process.poll()
            if pp != None:
                mylog.info("Process has died; stopping fido")
                return 1
        return 0


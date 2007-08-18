"""
Output-watcher

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
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
        if functionHandler == None:
            self.functionHandler = lambda a: None
        else:
            self.functionHandler = functionHandler()

    def run(self):
        mylog.info("Entering main Fido loop (CTRL-C or touch 'stopFido' to end)")
        while not self.checkForStop():
            nn = self.checkForOutput()
            for bn in nn:
                newName = buryOutput(bn)
                self.dealWithOutput(newName)
            time.sleep(WAITBETWEEN)

    def dealWithOutput(self, filename):
        # First, add it to the OutputCollection
        self.oc.addOutput(filename)
        self.oc.writeOut("runF_%s" % (self.title))
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
            filesFound = glob.glob("*.hierarchy")
            for file in filter(lambda a: a not in self.skipFiles, filesFound):
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


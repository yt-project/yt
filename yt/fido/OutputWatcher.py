"""
Output-watcher

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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
    def __init__(self, title=None, path=".", new_prefix="", oc=None,
                 process=None, function_handler = None):
        self.original_path = os.getcwd()
        os.chdir(path)
        if title == None: title = os.path.basename(os.getcwd())
        if oc == None: oc = OutputCollection(title)
        self.title = title
        self.oc = oc
        self.process = process
        self.new_prefix = new_prefix
        self.skip_files = [] # Forward compatible
        self.function_handler = None
        if function_handler != None:
            self.function_handler = function_handler()

    def run(self, run_once=False):
        mylog.info("Entering main Fido loop (CTRL-C or touch 'stopFido' to end)")
        wb = ytcfg.getfloat("fido","WaitBetween")
        while not self.check_for_stop():
            nn = self.check_for_output()
            for bn in nn:
                new_name = bury_output(bn, new_prefix=self.new_prefix)
                self.handle_output(new_name)
            if run_once: break
            try:
                time.sleep(wb)
            except KeyboardInterrupt:
                sys.exit()

    def handle_output(self, filename):
        # First, add it to the OutputCollection
        self.oc.add_output(filename)
        self.oc.write_out()
        if self.function_handler == None: return
        # Now, we pass it to our function handler
        pid = os.fork()
        if pid:
            mylog.debug("Waiting on pid %s", pid)
            newpid, exit = os.waitpid(pid,0)
            mylog.debug("Exit status %s from PID %s", exit, newpid)
        else:
            mylog.info("Forked process reporting for duty!")
            self.function_handler(filename)
            sys.exit()

    def check_for_output(self):
        new_files = []
        if os.path.isfile(NEW_OUTPUT_CREATED):
            os.unlink(NEW_OUTPUT_CREATED)
            # So something is created!  Now let's snag it
            # We insert our additional glob patterns here
            files_found = glob.glob("*.hierarchy") #\
                         #[glob.glob(pat) for pat in GlobPatterns]
            for file in [fn for fn in files_found if fn not in self.skip_files]:
                #if self.oc.has_key(os.path.basename(file[:-10])): continue
                new_files.append(file.rsplit(".",1)[0])
                mylog.info("Found output %s", new_files[-1])
        return new_files

    def check_for_stop(self):
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


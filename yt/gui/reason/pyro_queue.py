"""
Task queue to connect with reason via Pyro4.

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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

import numpy as na
import Pyro4
import uuid

from yt.funcs import *

from yt.gui.reason.basic_repl import ProgrammaticREPL
from yt.gui.reason.bottle_mods import PayloadHandler

class PyroQueueRoot(object):
    def __init__(self, comm):
        self.comm = comm
        self.repl = ProgrammaticREPL()
        self.payload_handler = PayloadHandler()

    def execute(self, code):
        mylog.info('Root sending out code.')
        code = self.comm.comm.bcast(code, root=0)
        value = self.repl.execute(code)
        return value

    def deliver(self):
        return self.payload_handler.deliver_payloads()

class PyroQueueNonRoot(object):
    def __init__(self, comm):
        self.comm = comm
        self.repl = ProgrammaticREPL()

    def run(self):
        while 1:
            code = None
            code = self.comm.comm.bcast(code, root=0)
            mylog.info('Received commands from subcomm root.')
            value = self.repl.execute(code)
            

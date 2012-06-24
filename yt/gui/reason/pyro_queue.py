"""
Task queue for reason.

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

from yt.funcs import *

from yt.gui.reason.basic_repl import ProgrammaticREPL
from yt.gui.reason.extdirect_repl import ExecutionThread
from yt.gui.reason.bottle_mods import PayloadHandler
from .utils import get_list_of_datasets

class PyroQueueRoot(object):
    def __init__(self, comm):
        self.comm = comm
        self.repl = ProgrammaticREPL()
        self.execution_thread = ExecutionThread(self.repl)
        self.payload_handler = PayloadHandler()
        self.execution_thread.start()

    def execute(self, code):
        mylog.info('Root sending out code.')
        code = self.comm.comm.bcast(code, root=0)
        self.execution_thread.execute_one(code, False)

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
            

"""
Task queue for reason.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

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

    def execute(self, code, hide = False):
        mylog.info('Root sending out code.')
        code = self.comm.comm.bcast(code, root=0)
        task = {'type': 'code',
                'code': code,
                'hide': hide}
        self.execution_thread.queue.put(task)

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
            

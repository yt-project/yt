"""
Task queue in yt

Author: Britton Smith <matthewturk@gmail.com>
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

from mpi4py import MPI
import time, threading, random

from .parallel_analysis_interface import \
    _get_comm

messages = dict(
    task = dict(msg = 'next'),
    result = dict(msg = 'result'),
    task_req = dict(msg = 'task_req'),
    end = dict(msg = 'no_more_tasks'),
)

class TaskQueueNonRoot(object):
    def __init__(self, tasks):
        self.tasks = tasks
        self.results = None
        self.comm = _get_comm(())

    def send_result(self, result):
        new_msg = messages['result'].copy()
        new_msg['value'] = result
        self.comm.comm.send(new_msg, dest = 0, tag=1)

    def get_next(self):
        msg = messages['task_req'].copy()
        self.comm.comm.send(msg, dest = 0, tag=1)
        msg = self.comm.comm.recv(source = 0, tag=2)
        if msg['msg'] == messages['end']['msg']:
            raise StopIteration
        return msg['value']

    def __iter__(self):
        while 1:
            yield self.get_next()

    def run(self, callable):
        for task in self:
            self.send_result(callable(task))
        return self.finalize()

    def finalize(self):
        return self.comm.comm.bcast(None, root = 0)

class TaskQueueRoot(TaskQueueNonRoot):
    def __init__(self, tasks):
        self.tasks = tasks
        self.results = {}
        self.assignments = {}
        self._notified = 0
        self._current = 0
        self._remaining = len(self.tasks)
        self.comm = _get_comm(())
        self.handle_assignments()
        # Set up threading here
        # self.dist = threading.Thread(target=self.handle_assignments)
        # self.dist.daemon = True
        # self.dist.start()

    def insert_result(self, source_id, result):
        task_id = self.assignments[source_id]
        self.results[task_id] = result

    def assign_task(self, source_id):
        if self._remaining == 0:
            msg = messages['end'].copy()
            self._notified += 1
        else:
            msg = messages['task'].copy()
            task_id = self._current
            task = self.tasks[task_id]
            self.assignments[source_id] = task_id
            self._current += 1
            self._remaining -= 1
            msg['value'] = task
        self.comm.comm.send(msg, dest = source_id, tag = 2)

    def handle_assignments(self):
        while 1:
            self.comm.probe(1, self.handshake)
            if self._notified >= self.comm.comm.size: break

    def handshake(self, status)
        msg = self.comm.recv(source = status.source, tag = 1)
        if msg['msg'] == messages['result']['msg']:
            self.insert_result(st.source, msg['value'])
        elif msg['msg'] == messages['task_req']['msg']:
            self.assign_task(st.source)
        else:
            print "GOT AN UNKNOWN MESSAGE", msg
            raise RuntimeError

    def finalize(self):
        self.dist.join()
        rv = self.comm.comm.bcast(self.results, root = 0)
        return rv


"""
A task queue for distributing work to worker agents

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

import threading
from yt.funcs import *

# The idea here is that we have a set of tasks, which we want to distribute.
# We'll try to make this forward-compatible.  To do so, we want to support the
# idea that there's a single, global set of tasks, as well as consumers that
# receive tasks from the main controller.  These consumers then pass them out
# to executors.
#
# The middle level, the "Consumer," is only really distinct from the executor
# in the case that there is an MPI subcommunicator.  The reason for the
# separation is so that the controller only communicates with a single member
# of each subcommunicator, which then passes that info back out.

def locked(func):
    @wraps(func)
    def exclusive(self, *args, **kwargs):
        with self.lock:
            return func(self, *args, **kwargs)
    return exclusive

class YTTaskCommunicator(object):
    # This should carefully be checked for a race condition, particularly in
    # the wait() function
    def __init__(self, interval = 2.0):
        self.interval = interval
        self.task_id = None
        self.waiting = False
        self.lock = threading.Lock()

    @locked
    def send_task(self, task_id):
        self.task_id = task_id

    @locked
    def query(self):
        return (self.waiting and self.task_id is None)

    def wait(self):
        self.waiting = True
        while self.task_id is None:
            time.sleep(self.interval)
        with self.lock:
            self.waiting = False
            new_task_id = self.task_id
            self.task_id = None
        return new_task_id

class YTTaskQueueController(threading.Thread):
    # There's only one of these for every instance of yt -- whether than
    # instance be spread across processors or not.
    # We assume that this will exist in the process space of a consumer, so it
    # will be threading based.
    def __init__(self, tasks, interval = 2.0, communicators = None):
        self.assignments = []
        self.interval = interval
        # Communicators can be anything but they have to implement a mechanism
        # for saying, "I'm ready" and "I'm done"
        self.tasks = tasks
        self.communicators = communicators
        threading.Thread.__init__(self)

    def run(self):
        # Now we bootstrap
        for i,c in enumerate(self.communicators):
            self.assignments.append(i)
            if i == len(self.tasks): break
            c.send_task(i)
        while len(self.assignments) < len(self.tasks):
            time.sleep(self.interval)
            for i,c in enumerate(self.communicators):
                if not c.query(): continue
                print "Sending assignment %s to %s" % (
                    len(self.assignments), i)
                c.send_task(len(self.assignments))
                self.assignments.append(i)
                if len(self.assignments) >= len(self.tasks): break
        terminated = 0
        while terminated != len(self.communicators):
            for i,c in enumerate(self.communicators):
                if not c.query(): continue
                c.send_task(-1)
                terminated += 1
                print "Terminated %s" % (i)

class YTTaskQueueConsumer(object):
    # One of these will exist per individual MPI task or one per MPI
    # subcommunicator, depending on the level of parallelism.  They serve to
    # talk to the YTTaskQueueController on one side and possibly several
    # YTTaskExecutors on the other.
    #
    # One potential setup for this, when using MPI, would be to have the
    # Executors each have one of these, but only the head process of that
    # subcommunicator possess an external communicator.  Then in next_task,
    # if the external communicator exists, one would probe that; otherwise,
    # accept a broadcast from the internal communicator's 0th task.
    def __init__(self, external_communicator, internal_communicator):
        self.external_communicator = external_communicator
        self.internal_communicator = internal_communicator

    def next_task(self):
        next_task = self.external_communicator.wait()
        #self.internal_communicator.notify(next_task)
        return next_task

class YTTaskExecutor(object):
    _count = 0
    # One of these will exist per computational actor
    def __init__(self, tasks, communicator):
        self.communicator = communicator
        self.tasks = tasks
        self.name = "Runner%03s" % (self.__class__._count)
        self.__class__._count += 1

    def run(self):
        # Note that right now this only works for a 1:1 mapping of
        # YTTaskQueueConsumer to YTTaskExecutor
        next_task = None
        while 1:
            next_task = self.communicator.next_task()
            if next_task == -1: break
            print "Executing on %s" % (self.name),
            self.tasks[next_task]()
        print "Concluded on %s" % (self.name)

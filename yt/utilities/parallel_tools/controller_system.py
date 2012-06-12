"""
A queueing system based on MPI

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia
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
    
try:
    from .parallel_analysis_interface import MPI
except ImportError:
    pass
from contextmanager import contextlib
from abc import ABCMeta, abstractmethod, abstractproperty

class WorkSplitter(object):
    def __init__(self, controller, group1, group2):
        self.group1 = group1
        self.group2 = group2
        self.controller = controller

    @classmethod
    def setup(cls, ng1, ng2):
        pp, wg = ProcessorPool.from_sizes(
            [(1, "controller"), (ng1, "group1"), (ng2, "group2")])
        groupc = pp['controller']
        group1 = pp['group1']
        group2 = pp['group2']
        obj = cls(groupc, group1, group2)
        obj.run(wg.name)

    def run(self, name):
        if name == "controller":
            self.run_controller()
        elif name == "group1":
            self.run_group1()
        elif name == "group2":
            self.run_group2()
        else:
            raise NotImplementedError

    @abstractmethod
    def run_controller(self):
        pass

    @abstractmethod
    def run_group1(self):
        pass

    @abstractmethod
    def run_group2(self):
        pass

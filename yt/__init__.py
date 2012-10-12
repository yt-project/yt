"""
YT is a package written primarily in Python designed to make the task of
running Enzo easier.  It contains facilities for creating Enzo data (currently
in prototype form) as well as runnning Enzo simulations, simulating the actions
of Enzo on various existing data, and analyzing output from Enzo in a
wide-variety of methods.

An ever-growing collection of documentation is also available at
http://yt-project.org/doc/ . Additionally, there is a
project site at http://yt-project.org/ with recipes, a wiki, a variety of
ways of peering into the version control, and a bug-reporting system.

YT is divided into several packages.

frontends
---------

This is where interfaces to codes are created.  Within each subdirectory of
yt/frontends/ there must exist the following files, even if empty:

* data_structures.py, where subclasses of AMRGridPatch, StaticOutput and
  AMRHierarchy are defined.
* io.py, where a subclass of IOHandler is defined.
* misc.py, where any miscellaneous functions or classes are defined.
* definitions.py, where any definitions specific to the frontend are
  defined.  (i.e., header formats, etc.)

visualization
-------------

This is where all visualization modules are stored.  This includes plot
collections, the volume rendering interface, and pixelization frontends.

data_objects
------------

All objects that handle data, processed or unprocessed, not explicitly
defined as visualization are located in here.  This includes the base
classes for data regions, covering grids, time series, and so on.  This
also includes derived fields and derived quantities.

analysis_modules
----------------

This is where all mechanisms for processing data live.  This includes
things like clump finding, halo profiling, halo finding, and so on.  This
is something of a catchall, but it serves as a level of greater
abstraction that simply data selection and modification.

gui
---

This is where all GUI components go.  Typically this will be some small
tool used for one or two things, which contains a launching mechanism on
the command line.

utilities
---------

All broadly useful code that doesn't clearly fit in one of the other
categories goes here.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

__version__ = "2.5-dev"

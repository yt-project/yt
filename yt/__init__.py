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

* data_structures.py, where subclasses of AMRGridPatch, Dataset and
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



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = "3.0-dev"

def run_nose(verbose=False, run_answer_tests=False, answer_big_data=False):
    import nose, os, sys
    from yt.config import ytcfg
    nose_argv = sys.argv
    nose_argv += ['--exclude=answer_testing','--detailed-errors']
    if verbose:
        nose_argv.append('-v')
    if run_answer_tests:
        nose_argv.append('--with-answer-testing')
    if answer_big_data:
        nose_argv.append('--answer-big-data')
    log_suppress = ytcfg.getboolean("yt","suppressStreamLogging")
    ytcfg.set("yt","suppressStreamLogging", 'True')
    initial_dir = os.getcwd()
    yt_file = os.path.abspath(__file__)
    yt_dir = os.path.dirname(yt_file)
    os.chdir(yt_dir)
    try:
        nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
        ytcfg.set("yt","suppressStreamLogging", str(log_suppress))

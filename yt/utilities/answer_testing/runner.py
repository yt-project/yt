"""
Runner mechanism for answer testing

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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

import matplotlib
import os, shelve, cPickle, sys, imp, tempfile

from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"
from yt.funcs import *
from yt.utilities.command_line import YTCommand
from .xunit import Xunit

from output_tests import test_registry, MultipleOutputTest, \
                         RegressionTestException

def clear_registry():
    test_registry.clear()

class FileNotExistException(Exception):
    def __init__(self, filename):
        self.filename = filename

    def __repr__(self):
        return "FileNotExistException: %s" % (self.filename)


def registry_entries():
    return test_registry.keys()

class RegressionTestStorage(object):
    def __init__(self, results_id, path = "."):
        self.id = results_id
        if results_id == "":
            self._path = os.path.join(path, "results")
        else:
            self._path = os.path.join(path, "results_%s" % self.id)
        if not os.path.isdir(self._path): 
            only_on_root(os.mkdir, self._path)
        if os.path.isfile(self._path): raise RuntimeError

    def _fn(self, tn):
        return os.path.join(self._path, tn)

    @rootonly
    def __setitem__(self, test_name, result):
        # We have to close our shelf manually,
        # as the destructor does not necessarily do this.
        # Context managers would be more appropriate.
        f = open(self._fn(test_name), "wb")
        cPickle.dump(result, f, protocol=-1)
        f.close()

    def __getitem__(self, test_name):
        if not os.path.exists(self._fn(test_name)):
            raise FileNotExistException(self._fn(test_name))
        f = open(self._fn(test_name), "rb")
        tr = cPickle.load(f)
        f.close()
        return tr

class RegressionTestRunner(object):
    def __init__(self, results_id, compare_id = None,
                 results_path = ".", compare_results_path = ".",
                 io_log = "OutputLog", plot_tests = False):
        # This test runner assumes it has been launched with the current
        # working directory that of the test case itself.
        self.io_log = io_log
        self.id = results_id
        if compare_id is not None:
            self.old_results = RegressionTestStorage(
                                    compare_id, path=compare_results_path)
        else:
            self.old_results = None
        self.results = RegressionTestStorage(results_id, path=results_path)
        self.plot_list = {}
        self.passed_tests = {}
        self.plot_tests = plot_tests

    def run_all_tests(self):
        plot_list = []
        for i,name in enumerate(sorted(test_registry)):
            self.run_test(name)
        return self.passed_tests

    def run_test(self, name):
        # We'll also need to call the "compare" operation,
        # but for that we'll need a data store.
        test = test_registry[name]
        plot_list = []
        if test.output_type == 'single':
            mot = MultipleOutputTest(self.io_log)
            for i,fn in enumerate(mot):
                # This next line is to keep the shelve module
                # from happily gobbling the disk
                #if i > 5: break 
                test_instance = test(fn)
                test_instance.name = "%s_%s" % (
                    os.path.basename(fn), test_instance.name )
                self._run(test_instance)

        elif test.output_type == 'multiple':
            test_instance = test(self.io_log)
            self._run(test_instance)

    watcher = None
    def _run(self, test):
        if self.watcher is not None:
            self.watcher.start()
        print self.id, "Running", test.name,
        test.setup()
        test.run()
        if self.plot_tests:
            self.plot_list[test.name] = test.plot()
        self.results[test.name] = test.result
        success, msg = self._compare(test)
        if self.old_results is None:
            print "NO OLD RESULTS"
        else:
            if success == True: print "SUCCEEDED"
            else: print "FAILED", msg
        self.passed_tests[test.name] = success
        if self.watcher is not None:
            if success == True:
                self.watcher.addSuccess(test.name)
            else:
                self.watcher.addFailure(test.name, msg)

    def _compare(self, test):
        if self.old_results is None:
            return (True, "New Test")
        try:
            old_result = self.old_results[test.name]
        except FileNotExistException:
            return (False, sys.exc_info())
        try:
            test.compare(old_result)
        except RegressionTestException, exc:
            return (False, sys.exc_info())
        return (True, "Pass")

    def run_tests_from_file(self, filename):
        for line in open(filename):
            test_name = line.strip()
            if test_name not in test_registry:
                if test_name[0] != "#":
                    print "Test '%s' not recognized, skipping" % (test_name)
                continue
            print "Running '%s'" % (test_name)
            self.run_test(line.strip())

def _load_modules(test_modules):
    for fn in test_modules:
        if fn.endswith(".py"): fn = fn[:-3]
        print "Loading module %s" % (fn)
        mname = os.path.basename(fn)
        f, filename, desc = imp.find_module(mname, [os.path.dirname(fn)])
        project = imp.load_module(mname, f, filename, desc)

def _update_io_log(opts, kwargs):
    if opts.datasets is None or len(opts.datasets) == 0: return
    f = tempfile.NamedTemporaryFile()
    kwargs['io_log'] = f.name
    for d in opts.datasets:
        fn = os.path.expanduser(d)
        print "Registered dataset %s" % fn
        f.write("DATASET WRITTEN %s\n" % fn)
    f.flush()
    f.seek(0)
    return f

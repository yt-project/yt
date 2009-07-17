"""
This is a definition of a class for defining problem verification.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from yt.mods import *
from yt.funcs import *
from yt.extensions.EnzoSimulation import *
from SimulationTests import test_registry

from mercurial import commands, hg, ui
import cPickle

class Problem(object):
    def __init__(self, name, pfname, repo_path="."):
        self.name = name
        self.pfname = pfname
        self.simulation = EnzoSimulation(pfname)
        self.tests, self.test_types = [], {}
        self.init_repo(repo_path)
        self._initialize_tests()

    def _add_existing_test(self, test_type):
        @wraps(test_type.__init__)
        def test_wrapper(name, *args, **kwargs):
            new_test = test_type(name, *args, **kwargs)
            self.tests.append(new_test)
        return test_wrapper

    def _initialize_tests(self):
        for name, test in test_registry.iteritems():
            setattr(self, 'add_%s' % name, self._add_existing_test(test))
            self.test_types[name] = self._add_existing_test(test)

    def init_repo(self, path):
        self.ui = ui.ui()
        try:
            self.repo = hg.repository(self.ui, path)
        except:
            print "CREATING:", path
            self.repo = hg.repository(self.ui, path, create=True)

    def load_repo_file(self, fn, identifier='tip'):
        return self.repo[identifier][fn].data()

    def run_tests(self):
        all_results = {}
        try:
            for pf in self.simulation:
                results = {}
                for test in self.tests:
                    results[test.name] = test(pf)
                all_results[str(pf)] = results
        except IOError:
            pass
        return all_results

    def store_results(self, results):
        base_fn = "results_%s.cpkl" % self.name
        fn = self.repo.pathto(base_fn)
        cPickle.dump(results, open(fn, "w"))
        if base_fn not in self.repo['tip'].manifest():
            commands.add(self.ui, self.repo, fn)
        message = "Committing results from current run"
        commands.commit(self.ui, self.repo, fn, message=message)
        print "Committed"

    def __call__(self):
        results = self.run_tests()
        self.store_results((self.tests, results))

# repo['tip']['yt/lagos/__init__.py'].data()

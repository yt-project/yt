"""
Answer Testing using Nose as a starting point

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

import logging
import os
import shelve
import hashlib
import contextlib

from yt.testing import *
from yt.utilities.command_line import get_yt_version
from yt.config import ytcfg
from nose.plugins import Plugin
from yt.mods import *

log = logging.getLogger('nose.plugins.answer-testing')

class AnswerTesting(Plugin):
    name = "answer-testing"

    def options(self, parser, env=os.environ):
        super(AnswerTesting, self).options(parser, env=env)
        test_storage_directory = ytcfg.get("yt", "test_storage_dir")
        try:
            my_hash = get_yt_version()
        except:
            my_hash = "UNKNOWN%s" % (time.time())
        parser.add_option("--answer-parameter-file", dest="parameter_file",
            default=os.path.join(os.getcwd(), "tests/DD0010/moving7_0010"),
            help="The parameter file value to feed to 'load' to test against")
        parser.add_option("--answer-output", dest="storage_dir",
            default=test_storage_directory,
            help="Base directory for storing test output.")
        parser.add_option("--answer-compare", dest="compare_name",
            default=None,
            help="The name against which we will compare")
        parser.add_option("--answer-name", dest="this_name",
            default=my_hash,
            help="The name we'll call this set of tests")

    def configure(self, options, conf):
        super(AnswerTesting, self).configure(options, conf)
        if not self.enabled:
            return
        AnswerTestingTest.result_storage = shelve.open(
            os.path.join(options.storage_dir,
                         options.this_name))
        if options.compare_name is not None:
            AnswerTestingTest.reference_storage = shelve.open(
                os.path.join(options.storage_dir,
                            options.compare_name))

    def finalize(self, result):
        pass

@contextlib.contextmanager
def temp_cwd(cwd):
    oldcwd = os.getcwd()
    os.chdir(cwd)
    yield
    os.chdir(oldcwd)

class AnswerTestingTest(object):
    reference_storage = None

    description = None
    def __init__(self, name, pf_fn):
        path = ytcfg.get("yt", "data_storage_dir")
        with temp_cwd(path):
            self.pf = load(pf_fn)
            self.pf.h
        self.name = "%s_%s" % (self.pf, name)

    def __call__(self):
        nv = self.run()
        self.result_storage[self.name] = nv
        if self.reference_storage is not None:
            ov = self.reference_storage.get(self.name, None)
            return self.compare(nv, ov)
        else:
            ov = None
            return True

    def compare(self, new_result, old_result):
        raise RuntimeError

    def create_obj(self, pf, obj_type):
        # obj_type should be tuple of
        #  ( obj_name, ( args ) )
        if obj_type is None:
            return pf.h.all_data()
        cls = getattr(pf.h, obj_type[0])
        obj = cls(*obj_type[1])
        return obj

    @property
    def sim_center(self):
        """
        This returns the center of the domain.
        """
        return 0.5*(self.pf.domain_right_edge + self.pf.domain_left_edge)

    @property
    def max_dens_location(self):
        """
        This is a helper function to return the location of the most dense
        point.
        """
        return self.pf.h.find_max("Density")[1]

    @property
    def entire_simulation(self):
        """
        Return an unsorted array of values that cover the entire domain.
        """
        return self.pf.h.all_data()
        
class FieldValuesTest(AnswerTestingTest):

    def __init__(self, name, pf_fn, field, obj_type = None):
        super(FieldValuesTest, self).__init__(name, pf_fn)
        self.obj_type = obj_type
        self.field = field

    def run(self):
        obj = self.create_obj(self.pf, self.obj_type)
        return obj[self.field].sort()

    def compare(self, new_result, old_result):
        assert_equal(new_result, old_result)

def _try_load(pf_fn):
    try:
        load(pf_fn)
    except:
        return False
    return True

def assert_fields(pf_fn, fields, data_obj = None):
    if AnswerTestingTest.result_storage is None: return 
    for field in fields:
        yield FieldValuesTest("FieldValues_%s" % field, pf_fn,
                              field, data_obj)

class ProjectionValuesTest(AnswerTestingTest):
    def __init__(self, name, pf_fn, axis, field,
                 weight_field = None, data_source = None):
        super(ProjectionValuesTest, self).__init__(name, pf_fn)
        self.axis = axis
        self.field = field
        self.weight_field = field
        self.data_source = None

    def run(self):
        proj = self.pf.h.proj(self.axis, self.field,
                              weight_field=self.weight_field)
        return proj

    def compare(self, new_result, old_result):
        assert(len(new_result.field_data) == len(old_result.field_data))
        for k in new_result.field_data:
            assert (k in old_result.field_data)
        for k in new_result:
            assert_equal(new_result[k], old_result[k])

class GridValuesTest(AnswerTestingTest):
    def __init__(self, name, pf_fn, field):
        super(GridValuesTest, self).__init__(name, pf_fn)

    def run(self):
        hashes = {}
        for g in self.pf.h.grids:
            hashes[g.id] = hashlib.md5(g[self.field].tostring()).hexdigest()
            g.clear_data()
        return hashes

    def compare(self, new_result, old_result):
        assert(len(new_result) == len(old_result))
        for k in new_result:
            assert (k in old_result)
        for k in new_result:
            assert_equal(new_result[k], old_result[k])

class TestGridHierarchy(AnswerTestingTest):
    def run(self):
        result = {}
        result["grid_dimensions"] = self.pf.h.grid_dimensions
        result["grid_left_edges"] = self.pf.h.grid_left_edge
        result["grid_right_edges"] = self.pf.h.grid_right_edge
        result["grid_levels"] = self.pf.h.grid_levels
        result["grid_particle_count"] = self.pf.h.grid_particle_count

    def compare(self, new_result, old_result):
        for k in new_result:
            assert_equal(new_result[k], old_result[k])

class TestParentageRelationships(AnswerTestingTest):
    def run(self):
        result = {}
        result["parents"] = []
        result["children"] = []
        for g in self.pf.h.grids:
            p = g.Parent
            if p is None:
                result["parents"].append(None)
            elif hasattr(p, "id"):
                result["parents"].append(p.id)
            else:
                result["parents"].append([pg.id for pg in p])
            result["children"].append([c.id for c in g.Children])
        return result

    def compare(self, new_result, old_result):
        for newp, oldp in zip(new_result["parents"], old_result["parents"]):
            assert(newp == oldp)
        for newc, oldc in zip(new_result["children"], old_result["children"]):
            assert(newp == oldp)

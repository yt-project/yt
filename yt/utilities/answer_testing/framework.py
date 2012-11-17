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
import hashlib
import contextlib
import urllib2
import cPickle
import sys

from nose.plugins import Plugin
from yt.testing import *
from yt.config import ytcfg
from yt.mods import *
from yt.data_objects.static_output import StaticOutput
import cPickle
import shelve

from yt.utilities.logger import disable_stream_logging
from yt.utilities.command_line import get_yt_version

mylog = logging.getLogger('nose.plugins.answer-testing')
run_big_data = False

_latest = "gold001"
_url_path = "http://yt-answer-tests.s3-website-us-east-1.amazonaws.com/%s_%s"

class AnswerTesting(Plugin):
    name = "answer-testing"
    _my_version = None

    def options(self, parser, env=os.environ):
        super(AnswerTesting, self).options(parser, env=env)
        parser.add_option("--answer-compare-name", dest="compare_name", metavar='str',
            default=_latest, help="The name of tests against which we will compare")
        parser.add_option("--answer-big-data", dest="big_data",
            default=False, help="Should we run against big data, too?",
            action="store_true")
        parser.add_option("--answer-store-name", dest="store_name", metavar='str',
            default=None,
            help="The name we'll call this set of tests")
        parser.add_option("--local-store", dest="store_local_results",
            default=False, action="store_true", help="Store/Load local results?")

    @property
    def my_version(self, version=None):
        if self._my_version is not None:
            return self._my_version
        if version is None:
            try:
                version = get_yt_version()
            except:
                version = "UNKNOWN%s" % (time.time())
        self._my_version = version
        return self._my_version

    def configure(self, options, conf):
        super(AnswerTesting, self).configure(options, conf)
        if not self.enabled:
            return
        disable_stream_logging()
        if options.store_name is not None:
            self.store_results = True
        # If the user sets the storage_name, then it means they are storing and
        # not comparing, even if they set the compare_name (since it is set by default)
            options.compare_name = None
        else: 
            self.store_results = False
            options.store_name = self.my_version
        from yt.config import ytcfg
        ytcfg["yt","__withintesting"] = "True"
        AnswerTestingTest.result_storage = \
            self.result_storage = defaultdict(dict)
        if options.compare_name == "SKIP":
            options.compare_name = None
        elif options.compare_name == "latest":
            options.compare_name = _latest
            
        # Local/Cloud storage 
        if options.store_local_results:
            storage_class = AnswerTestLocalStorage
            # Fix up filename for local storage 
            if options.compare_name is not None:
                options.compare_name = "%s/%s/%s" % \
                    (os.path.realpath(options.output_dir), options.compare_name, 
                     options.compare_name)
            if options.store_name is not None:
                name_dir_path = "%s/%s" % \
                    (os.path.realpath(options.output_dir), 
                    options.store_name)
                if not os.path.isdir(name_dir_path):
                    os.makedirs(name_dir_path)
                options.store_name= "%s/%s" % \
                        (name_dir_path, options.store_name)
        else:
            storage_class = AnswerTestCloudStorage

        # Initialize answer/reference storage
        AnswerTestingTest.reference_storage = self.storage = \
                storage_class(options.compare_name, options.store_name)

        self.store_local_results = options.store_local_results
        global run_big_data
        run_big_data = options.big_data

    def finalize(self, result=None):
        if self.store_results is False: return
        self.storage.dump(self.result_storage)        

class AnswerTestStorage(object):
    def __init__(self, reference_name=None, answer_name=None):
        self.reference_name = reference_name
        self.answer_name = answer_name
        self.cache = {}
    def dump(self, result_storage, result):
        raise NotImplementedError 
    def get(self, pf_name, default=None):
        raise NotImplementedError 

class AnswerTestCloudStorage(AnswerTestStorage):
    def get(self, pf_name, default = None):
        if self.reference_name is None: return default
        if pf_name in self.cache: return self.cache[pf_name]
        url = _url_path % (self.reference_name, pf_name)
        try:
            resp = urllib2.urlopen(url)
            # This is dangerous, but we have a controlled S3 environment
            data = resp.read()
            rv = cPickle.loads(data)
        except urllib2.HTTPError as ex:
            raise YTNoOldAnswer(url)
            mylog.warning("Missing %s (%s)", url, ex)
            rv = default
        self.cache[pf_name] = rv
        return rv

    def dump(self, result_storage):
        if self.answer_name is None: return
        # This is where we dump our result storage up to Amazon, if we are able
        # to.
        import boto
        from boto.s3.key import Key
        c = boto.connect_s3()
        bucket = c.get_bucket("yt-answer-tests")
        for pf_name in result_storage:
            rs = cPickle.dumps(result_storage[pf_name])
            tk = bucket.get_key("%s_%s" % (self.answer_name, pf_name)) 
            if tk is not None: tk.delete()
            k = Key(bucket)
            k.key = "%s_%s" % (self.answer_name, pf_name)
            k.set_contents_from_string(rs)
            k.set_acl("public-read")

class AnswerTestLocalStorage(AnswerTestStorage):
    def dump(self, result_storage):
        if self.answer_name is None: return
        # Store data using shelve
        ds = shelve.open(self.answer_name, protocol=-1)
        for pf_name in result_storage:
            answer_name = "%s" % pf_name
            if name in ds:
                mylog.info("Overwriting %s", answer_name)
            ds[answer_name] = result_storage[pf_name]
        ds.close()

    def get(self, pf_name, default=None):
        if self.reference_name is None: return default
        # Read data using shelve
        answer_name = "%s" % pf_name
        ds = shelve.open(self.reference_name, protocol=-1)
        try:
            result = ds[answer_name]
        except KeyError:
            result = default
        ds.close()
        return result

@contextlib.contextmanager
def temp_cwd(cwd):
    oldcwd = os.getcwd()
    os.chdir(cwd)
    yield
    os.chdir(oldcwd)

def can_run_pf(pf_fn):
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    if isinstance(pf_fn, StaticOutput):
        return AnswerTestingTest.result_storage is not None
    with temp_cwd(path):
        try:
            load(pf_fn)
        except:
            return False
    return AnswerTestingTest.result_storage is not None

def data_dir_load(pf_fn):
    path = ytcfg.get("yt", "test_data_dir")
    if isinstance(pf_fn, StaticOutput): return pf_fn
    if not os.path.isdir(path):
        return False
    with temp_cwd(path):
        pf = load(pf_fn)
        pf.h
        return pf

def sim_dir_load(sim_fn, path = None, sim_type = "Enzo",
                 find_outputs=False):
    if path is None and not os.path.exists(sim_fn):
        raise IOError
    if os.path.exists(sim_fn) or not path:
        path = "."
    with temp_cwd(path):
        return simulation(sim_fn, sim_type,
                          find_outputs=find_outputs)

class AnswerTestingTest(object):
    reference_storage = None
    result_storage = None
    prefix = ""
    def __init__(self, pf_fn):
        self.pf = data_dir_load(pf_fn)

    def __call__(self):
        nv = self.run()
        if self.reference_storage.reference_name is not None:
            dd = self.reference_storage.get(self.storage_name)
            if dd is None: raise YTNoOldAnswer(self.storage_name)
            ov = dd[self.description]
            self.compare(nv, ov)
        else:
            ov = None
        self.result_storage[self.storage_name][self.description] = nv

    @property
    def storage_name(self):
        if self.prefix != "":
            return "%s_%s" % (self.prefix, self.pf)
        return str(self.pf)

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

    @property
    def description(self):
        obj_type = getattr(self, "obj_type", None)
        if obj_type is None:
            oname = "all"
        else:
            oname = "_".join((str(s) for s in obj_type))
        args = [self._type_name, str(self.pf), oname]
        args += [str(getattr(self, an)) for an in self._attrs]
        return "_".join(args)
        
class FieldValuesTest(AnswerTestingTest):
    _type_name = "FieldValues"
    _attrs = ("field", )

    def __init__(self, pf_fn, field, obj_type = None,
                 decimals = None):
        super(FieldValuesTest, self).__init__(pf_fn)
        self.obj_type = obj_type
        self.field = field
        self.decimals = decimals

    def run(self):
        obj = self.create_obj(self.pf, self.obj_type)
        avg = obj.quantities["WeightedAverageQuantity"](self.field,
                             weight="Ones")
        (mi, ma), = obj.quantities["Extrema"](self.field)
        return np.array([avg, mi, ma])

    def compare(self, new_result, old_result):
        err_msg = "Field values for %s not equal." % self.field
        if self.decimals is None:
            assert_equal(new_result, old_result, 
                         err_msg=err_msg, verbose=True)
        else:
            assert_allclose(new_result, old_result, 10.**(-self.decimals),
                             err_msg=err_msg, verbose=True)

class AllFieldValuesTest(AnswerTestingTest):
    _type_name = "AllFieldValues"
    _attrs = ("field", )

    def __init__(self, pf_fn, field, obj_type = None,
                 decimals = None):
        super(AllFieldValuesTest, self).__init__(pf_fn)
        self.obj_type = obj_type
        self.field = field
        self.decimals = decimals

    def run(self):
        obj = self.create_obj(self.pf, self.obj_type)
        return obj[self.field]

    def compare(self, new_result, old_result):
        err_msg = "All field values for %s not equal." % self.field
        if self.decimals is None:
            assert_equal(new_result, old_result, 
                         err_msg=err_msg, verbose=True)
        else:
            assert_rel_equal(new_result, old_result, self.decimals,
                             err_msg=err_msg, verbose=True)
            
class ProjectionValuesTest(AnswerTestingTest):
    _type_name = "ProjectionValues"
    _attrs = ("field", "axis", "weight_field")

    def __init__(self, pf_fn, axis, field, weight_field = None,
                 obj_type = None, decimals = None):
        super(ProjectionValuesTest, self).__init__(pf_fn)
        self.axis = axis
        self.field = field
        self.weight_field = field
        self.obj_type = obj_type
        self.decimals = decimals

    def run(self):
        if self.obj_type is not None:
            obj = self.create_obj(self.pf, self.obj_type)
        else:
            obj = None
        proj = self.pf.h.proj(self.axis, self.field,
                              weight_field=self.weight_field,
                              data_source = obj)
        return proj.field_data

    def compare(self, new_result, old_result):
        assert(len(new_result) == len(old_result))
        for k in new_result:
            assert (k in old_result)
        for k in new_result:
            err_msg = "%s values of %s (%s weighted) projection (axis %s) not equal." % \
              (k, self.field, self.weight_field, self.axis)
            if self.decimals is None:
                assert_equal(new_result[k], old_result[k],
                             err_msg=err_msg)
            else:
                assert_allclose(new_result[k], old_result[k], 
                                 10.**-(self.decimals), err_msg=err_msg)

class PixelizedProjectionValuesTest(AnswerTestingTest):
    _type_name = "PixelizedProjectionValues"
    _attrs = ("field", "axis", "weight_field")

    def __init__(self, pf_fn, axis, field, weight_field = None,
                 obj_type = None):
        super(PixelizedProjectionValuesTest, self).__init__(pf_fn)
        self.axis = axis
        self.field = field
        self.weight_field = field
        self.obj_type = obj_type

    def run(self):
        if self.obj_type is not None:
            obj = self.create_obj(self.pf, self.obj_type)
        else:
            obj = None
        proj = self.pf.h.proj(self.axis, self.field,
                              weight_field=self.weight_field,
                              data_source = obj)
        frb = proj.to_frb((1.0, 'unitary'), 256)
        frb[self.field]
        frb[self.weight_field]
        d = frb.data
        d.update( dict( (("%s_sum" % f, proj[f].sum(dtype="float64"))
                         for f in proj.field_data.keys()) ) )
        return d

    def compare(self, new_result, old_result):
        assert(len(new_result) == len(old_result))
        for k in new_result:
            assert (k in old_result)
        for k in new_result:
            assert_rel_equal(new_result[k], old_result[k], 10)

class GridValuesTest(AnswerTestingTest):
    _type_name = "GridValues"
    _attrs = ("field",)

    def __init__(self, pf_fn, field):
        super(GridValuesTest, self).__init__(pf_fn)
        self.field = field

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

class VerifySimulationSameTest(AnswerTestingTest):
    _type_name = "VerifySimulationSame"
    _attrs = ()

    def __init__(self, simulation_obj):
        self.pf = simulation_obj

    def run(self):
        result = [ds.current_time for ds in self.pf]
        return result

    def compare(self, new_result, old_result):
        assert_equal(len(new_result), len(old_result),
                     err_msg="Number of outputs not equal.",
                     verbose=True)
        for i in range(len(new_result)):
            assert_equal(new_result[i], old_result[i],
                         err_msg="Output times not equal.",
                         verbose=True)
        
class GridHierarchyTest(AnswerTestingTest):
    _type_name = "GridHierarchy"
    _attrs = ()

    def run(self):
        result = {}
        result["grid_dimensions"] = self.pf.h.grid_dimensions
        result["grid_left_edges"] = self.pf.h.grid_left_edge
        result["grid_right_edges"] = self.pf.h.grid_right_edge
        result["grid_levels"] = self.pf.h.grid_levels
        result["grid_particle_count"] = self.pf.h.grid_particle_count
        return result

    def compare(self, new_result, old_result):
        for k in new_result:
            assert_equal(new_result[k], old_result[k])

class ParentageRelationshipsTest(AnswerTestingTest):
    _type_name = "ParentageRelationships"
    _attrs = ()
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

def requires_pf(pf_fn, big_data = False):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data == False and big_data == True:
        return ffalse
    elif not can_run_pf(pf_fn):
        return ffalse
    else:
        return ftrue

def small_patch_amr(pf_fn, fields):
    if not can_run_pf(pf_fn): return
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    yield GridHierarchyTest(pf_fn)
    yield ParentageRelationshipsTest(pf_fn)
    for field in fields:
        yield GridValuesTest(pf_fn, field)
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield ProjectionValuesTest(
                        pf_fn, axis, field, weight_field,
                        ds)
                yield FieldValuesTest(
                        pf_fn, field, ds)

def big_patch_amr(pf_fn, fields):
    if not can_run_pf(pf_fn): return
    dso = [ None, ("sphere", ("max", (0.1, 'unitary')))]
    yield GridHierarchyTest(pf_fn)
    yield ParentageRelationshipsTest(pf_fn)
    for field in fields:
        yield GridValuesTest(pf_fn, field)
        for axis in [0, 1, 2]:
            for ds in dso:
                for weight_field in [None, "Density"]:
                    yield PixelizedProjectionValuesTest(
                        pf_fn, axis, field, weight_field,
                        ds)

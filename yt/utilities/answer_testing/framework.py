"""
Answer Testing using Nose as a starting point



"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import logging
import numpy as np
import os
import time
import hashlib
import contextlib
import sys
from yt.extern.six.moves import cPickle, urllib
import shelve
import zlib
import tempfile
import glob

from collections import defaultdict

from matplotlib.testing.compare import compare_images
from nose.plugins import Plugin
from yt.funcs import \
    get_pbar
from yt.testing import \
    assert_equal, \
    assert_allclose_units, \
    assert_rel_equal, \
    assert_almost_equal
from yt.convenience import load, simulation
from yt.config import ytcfg
from yt.data_objects.static_output import Dataset
from yt.data_objects.time_series import SimulationTimeSeries
from yt.utilities.exceptions import \
    YTNoOldAnswer, \
    YTCloudError, \
    YTOutputNotIdentified, \
    YTNoAnswerNameSpecified
from yt.utilities.logger import disable_stream_logging
from yt.utilities.command_line import get_yt_version

import matplotlib.image as mpimg
import yt.visualization.plot_window as pw
import yt.visualization.particle_plots as particle_plots
import yt.visualization.profile_plotter as profile_plotter

mylog = logging.getLogger('nose.plugins.answer-testing')
run_big_data = False

# Set the latest gold and local standard filenames
_latest = ytcfg.get("yt", "gold_standard_filename")
_latest_local = ytcfg.get("yt", "local_standard_filename")
_url_path = ytcfg.get("yt", "answer_tests_url")

class AnswerTesting(Plugin):
    name = "answer-testing"
    _my_version = None

    def options(self, parser, env=os.environ):
        super(AnswerTesting, self).options(parser, env=env)
        parser.add_option("--answer-name", dest="answer_name", metavar='str',
            default=None, help="The name of the standard to store/compare against")
        parser.add_option("--answer-store", dest="store_results", metavar='bool',
            default=False, action="store_true",
            help="Should we store this result instead of comparing?")
        parser.add_option("--local", dest="local_results",
            default=False, action="store_true", help="Store/load reference results locally?")
        parser.add_option("--answer-big-data", dest="big_data",
            default=False, help="Should we run against big data, too?",
            action="store_true")
        parser.add_option("--local-dir", dest="output_dir", metavar='str',
                          help="The name of the directory to store local results")

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

        # Parse through the storage flags to make sense of them
        # and use reasonable defaults
        # If we're storing the data, default storage name is local
        # latest version
        if options.store_results:
            if options.answer_name is None:
                self.store_name = _latest_local
            else:
                self.store_name = options.answer_name
            self.compare_name = None
        # if we're not storing, then we're comparing, and we want default
        # comparison name to be the latest gold standard
        # either on network or local
        else:
            if options.answer_name is None:
                if options.local_results:
                    self.compare_name = _latest_local
                else:
                    self.compare_name = _latest
            else:
                self.compare_name = options.answer_name
            self.store_name = self.my_version

        self.store_results = options.store_results

        ytcfg["yt","__withintesting"] = "True"
        AnswerTestingTest.result_storage = \
            self.result_storage = defaultdict(dict)
        if self.compare_name == "SKIP":
            self.compare_name = None
        elif self.compare_name == "latest":
            self.compare_name = _latest

        # Local/Cloud storage
        if options.local_results:
            if options.output_dir is None:
                print('Please supply an output directory with the --local-dir option')
                sys.exit(1)
            storage_class = AnswerTestLocalStorage
            output_dir = os.path.realpath(options.output_dir)
            # Fix up filename for local storage
            if self.compare_name is not None:
                self.compare_name = os.path.join(output_dir, self.compare_name,
                                                 self.compare_name)

            # Create a local directory only when `options.answer_name` is
            # provided. If it is not provided then creating local directory
            # will depend on the `AnswerTestingTest.answer_name` value of the
            # test, this case is handled in AnswerTestingTest class.
            if options.store_results and options.answer_name is not None:
                name_dir_path = os.path.join(output_dir, self.store_name)
                if not os.path.isdir(name_dir_path):
                    os.makedirs(name_dir_path)
                self.store_name = os.path.join(name_dir_path, self.store_name)
        else:
            storage_class = AnswerTestCloudStorage

        # Initialize answer/reference storage
        AnswerTestingTest.reference_storage = self.storage = \
                storage_class(self.compare_name, self.store_name)
        AnswerTestingTest.options = options

        self.local_results = options.local_results
        global run_big_data
        run_big_data = options.big_data

    def finalize(self, result=None):
        if self.store_results is False: return
        self.storage.dump(self.result_storage)

    def help(self):
        return "yt answer testing support"

class AnswerTestStorage(object):
    def __init__(self, reference_name=None, answer_name=None):
        self.reference_name = reference_name
        self.answer_name = answer_name
        self.cache = {}
    def dump(self, result_storage, result):
        raise NotImplementedError
    def get(self, ds_name, default=None):
        raise NotImplementedError

class AnswerTestCloudStorage(AnswerTestStorage):
    def get(self, ds_name, default = None):
        if self.reference_name is None: return default
        if ds_name in self.cache: return self.cache[ds_name]
        url = _url_path.format(self.reference_name, ds_name)
        try:
            resp = urllib.request.urlopen(url)
        except urllib.error.HTTPError:
            raise YTNoOldAnswer(url)
        else:
            for this_try in range(3):
                try:
                    data = resp.read()
                except:
                    time.sleep(0.01)
                else:
                    # We were succesful
                    break
            else:
                # Raise error if all tries were unsuccessful
                raise YTCloudError(url)
            # This is dangerous, but we have a controlled S3 environment
            rv = cPickle.loads(data)
        self.cache[ds_name] = rv
        return rv

    def progress_callback(self, current, total):
        self.pbar.update(current)

    def dump(self, result_storage):
        if self.answer_name is None: return
        # This is where we dump our result storage up to Amazon, if we are able
        # to.
        import pyrax
        credentials = os.path.expanduser(os.path.join('~', '.yt', 'rackspace'))
        pyrax.set_credential_file(credentials)
        cf = pyrax.cloudfiles
        c = cf.get_container("yt-answer-tests")
        pb = get_pbar("Storing results ", len(result_storage))
        for i, ds_name in enumerate(result_storage):
            pb.update(i)
            rs = cPickle.dumps(result_storage[ds_name])
            object_name = "%s_%s" % (self.answer_name, ds_name)
            if object_name in c.get_object_names():
                obj = c.get_object(object_name)
                c.delete_object(obj)
            c.store_object(object_name, rs)
        pb.finish()

class AnswerTestLocalStorage(AnswerTestStorage):
    def dump(self, result_storage):
        # The 'tainted' attribute is automatically set to 'True'
        # if the dataset required for an answer test is missing
        # (see can_run_ds() and can_run_sim()).
        # This logic check prevents creating a shelve with empty answers.
        storage_is_tainted = result_storage.get('tainted', False)
        if self.answer_name is None or storage_is_tainted:
            return
        # Store data using shelve
        ds = shelve.open(self.answer_name, protocol=-1)
        for ds_name in result_storage:
            answer_name = "%s" % ds_name
            if answer_name in ds:
                mylog.info("Overwriting %s", answer_name)
            ds[answer_name] = result_storage[ds_name]
        ds.close()

    def get(self, ds_name, default=None):
        if self.reference_name is None: return default
        # Read data using shelve
        answer_name = "%s" % ds_name
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

def can_run_ds(ds_fn, file_check = False):
    result_storage = AnswerTestingTest.result_storage
    if isinstance(ds_fn, Dataset):
        return result_storage is not None
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    if file_check:
        return os.path.isfile(os.path.join(path, ds_fn)) and \
            result_storage is not None
    try:
        load(ds_fn)
    except YTOutputNotIdentified:
        if ytcfg.getboolean("yt", "requires_ds_strict"):
            if result_storage is not None:
                result_storage['tainted'] = True
            raise
        return False
    return result_storage is not None

def can_run_sim(sim_fn, sim_type, file_check = False):
    result_storage = AnswerTestingTest.result_storage
    if isinstance(sim_fn, SimulationTimeSeries):
        return result_storage is not None
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    if file_check:
        return os.path.isfile(os.path.join(path, sim_fn)) and \
            result_storage is not None
    try:
        simulation(sim_fn, sim_type)
    except YTOutputNotIdentified:
        if ytcfg.getboolean("yt", "requires_ds_strict"):
            if result_storage is not None:
                result_storage['tainted'] = True
            raise
        return False
    return result_storage is not None

def data_dir_load(ds_fn, cls = None, args = None, kwargs = None):
    args = args or ()
    kwargs = kwargs or {}
    path = ytcfg.get("yt", "test_data_dir")
    if isinstance(ds_fn, Dataset): return ds_fn
    if not os.path.isdir(path):
        return False
    if cls is None:
        ds = load(ds_fn, *args, **kwargs)
    else:
        ds = cls(os.path.join(path, ds_fn), *args, **kwargs)
    ds.index
    return ds

def sim_dir_load(sim_fn, path = None, sim_type = "Enzo",
                 find_outputs=False):
    if path is None and not os.path.exists(sim_fn):
        raise IOError
    if os.path.exists(sim_fn) or not path:
        path = "."
    return simulation(os.path.join(path, sim_fn), sim_type,
                      find_outputs=find_outputs)

class AnswerTestingTest(object):
    reference_storage = None
    result_storage = None
    prefix = ""
    options = None
    # This variable should be set if we are not providing `--answer-name` as
    # command line parameter while running yt's answer testing using nosetests.
    answer_name = None
    def __init__(self, ds_fn):
        if ds_fn is None:
            self.ds = None
        elif isinstance(ds_fn, Dataset):
            self.ds = ds_fn
        else:
            self.ds = data_dir_load(ds_fn)

    def __call__(self):
        if AnswerTestingTest.result_storage is None:
            return
        nv = self.run()

        # Test answer name should be provided either as command line parameters
        # or by setting AnswerTestingTest.answer_name
        if self.options.answer_name is None and self.answer_name is None:
            raise YTNoAnswerNameSpecified()

        # This is for running answer test when `--answer-name` is not set in
        # nosetests command line arguments. In this case, set the answer_name
        # from the `answer_name` keyword in the test case
        if self.options.answer_name is None:
            pyver = "py{}{}".format(sys.version_info.major,
                                    sys.version_info.minor)
            self.answer_name = "{}_{}".format(pyver, self.answer_name)

            answer_store_dir = os.path.realpath(self.options.output_dir)
            ref_name = os.path.join(answer_store_dir, self.answer_name,
                                    self.answer_name)
            self.reference_storage.reference_name = ref_name
            self.reference_storage.answer_name = ref_name

            # If we are generating golden answers (passed --answer-store arg):
            # - create the answer directory for this test
            # - self.reference_storage.answer_name will be path to answer files
            if self.options.store_results:
                answer_test_dir = os.path.join(answer_store_dir,
                                               self.answer_name)
                if not os.path.isdir(answer_test_dir):
                    os.makedirs(answer_test_dir)
                self.reference_storage.reference_name = None

        if self.reference_storage.reference_name is not None:
            # Compare test generated values against the golden answer
            dd = self.reference_storage.get(self.storage_name)
            if dd is None or self.description not in dd:
                raise YTNoOldAnswer(
                    "%s : %s" % (self.storage_name, self.description))
            ov = dd[self.description]
            self.compare(nv, ov)
        else:
            # Store results, hence do nothing (in case of --answer-store arg)
            ov = None
        self.result_storage[self.storage_name][self.description] = nv

    @property
    def storage_name(self):
        if self.prefix != "":
            return "%s_%s" % (self.prefix, self.ds)
        return str(self.ds)

    def compare(self, new_result, old_result):
        raise RuntimeError

    def create_plot(self, ds, plot_type, plot_field, plot_axis, plot_kwargs = None):
        # plot_type should be a string
        # plot_kwargs should be a dict
        if plot_type is None:
            raise RuntimeError('Must explicitly request a plot type')
        cls = getattr(pw, plot_type, None)
        if cls is None:
            cls = getattr(particle_plots, plot_type)
        plot = cls(*(ds, plot_axis, plot_field), **plot_kwargs)
        return plot

    @property
    def sim_center(self):
        """
        This returns the center of the domain.
        """
        return 0.5*(self.ds.domain_right_edge + self.ds.domain_left_edge)

    @property
    def max_dens_location(self):
        """
        This is a helper function to return the location of the most dense
        point.
        """
        return self.ds.find_max("density")[1]

    @property
    def entire_simulation(self):
        """
        Return an unsorted array of values that cover the entire domain.
        """
        return self.ds.all_data()

    @property
    def description(self):
        obj_type = getattr(self, "obj_type", None)
        if obj_type is None:
            oname = "all"
        else:
            oname = "_".join((str(s) for s in obj_type))
        args = [self._type_name, str(self.ds), oname]
        args += [str(getattr(self, an)) for an in self._attrs]
        return "_".join(args).replace('.', '_')

class FieldValuesTest(AnswerTestingTest):
    _type_name = "FieldValues"
    _attrs = ("field", )

    def __init__(self, ds_fn, field, obj_type = None,
                 particle_type=False, decimals = 10):
        super(FieldValuesTest, self).__init__(ds_fn)
        self.obj_type = obj_type
        self.field = field
        self.particle_type = particle_type
        self.decimals = decimals

    def run(self):
        obj = create_obj(self.ds, self.obj_type)
        field = obj._determine_fields(self.field)[0]
        if self.particle_type:
            weight_field = (field[0], "particle_ones")
        else:
            weight_field = ("index", "ones")
        avg = obj.quantities.weighted_average_quantity(
            field, weight=weight_field)
        mi, ma = obj.quantities.extrema(self.field)
        return np.array([avg, mi, ma])

    def compare(self, new_result, old_result):
        err_msg = "Field values for %s not equal." % (self.field,)
        if self.decimals is None:
            assert_equal(new_result, old_result,
                         err_msg=err_msg, verbose=True)
        else:
            assert_allclose_units(new_result, old_result, 10.**(-self.decimals),
                                  err_msg=err_msg, verbose=True)

class AllFieldValuesTest(AnswerTestingTest):
    _type_name = "AllFieldValues"
    _attrs = ("field", )

    def __init__(self, ds_fn, field, obj_type = None,
                 decimals = None):
        super(AllFieldValuesTest, self).__init__(ds_fn)
        self.obj_type = obj_type
        self.field = field
        self.decimals = decimals

    def run(self):
        obj = create_obj(self.ds, self.obj_type)
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

    def __init__(self, ds_fn, axis, field, weight_field = None,
                 obj_type = None, decimals = None):
        super(ProjectionValuesTest, self).__init__(ds_fn)
        self.axis = axis
        self.field = field
        self.weight_field = weight_field
        self.obj_type = obj_type
        self.decimals = decimals

    def run(self):
        if self.obj_type is not None:
            obj = create_obj(self.ds, self.obj_type)
        else:
            obj = None
        if self.ds.domain_dimensions[self.axis] == 1: return None
        proj = self.ds.proj(self.field, self.axis,
                              weight_field=self.weight_field,
                              data_source = obj)
        return proj.field_data

    def compare(self, new_result, old_result):
        if new_result is None:
            return
        assert(len(new_result) == len(old_result))
        nind, oind = None, None
        for k in new_result:
            assert (k in old_result)
            if oind is None:
                oind = np.array(np.isnan(old_result[k]))
            np.logical_or(oind, np.isnan(old_result[k]), oind)
            if nind is None:
                nind = np.array(np.isnan(new_result[k]))
            np.logical_or(nind, np.isnan(new_result[k]), nind)
        oind = ~oind
        nind = ~nind
        for k in new_result:
            err_msg = "%s values of %s (%s weighted) projection (axis %s) not equal." % \
              (k, self.field, self.weight_field, self.axis)
            if k == 'weight_field' and self.weight_field is None:
                continue
            nres, ores = new_result[k][nind], old_result[k][oind]
            if self.decimals is None:
                assert_equal(nres, ores, err_msg=err_msg)
            else:
                assert_allclose_units(nres, ores, 10.**-(self.decimals),
                                      err_msg=err_msg)

class PixelizedProjectionValuesTest(AnswerTestingTest):
    _type_name = "PixelizedProjectionValues"
    _attrs = ("field", "axis", "weight_field")

    def __init__(self, ds_fn, axis, field, weight_field = None,
                 obj_type = None):
        super(PixelizedProjectionValuesTest, self).__init__(ds_fn)
        self.axis = axis
        self.field = field
        self.weight_field = weight_field
        self.obj_type = obj_type

    def run(self):
        if self.obj_type is not None:
            obj = create_obj(self.ds, self.obj_type)
        else:
            obj = None
        proj = self.ds.proj(self.field, self.axis,
                              weight_field=self.weight_field,
                              data_source = obj)
        frb = proj.to_frb((1.0, 'unitary'), 256)
        frb[self.field]
        if self.weight_field is not None:
            frb[self.weight_field]
        d = frb.data
        for f in proj.field_data:
            # Sometimes f will be a tuple.
            d["%s_sum" % (f,)] = proj.field_data[f].sum(dtype="float64")
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

    def __init__(self, ds_fn, field):
        super(GridValuesTest, self).__init__(ds_fn)
        self.field = field

    def run(self):
        hashes = {}
        for g in self.ds.index.grids:
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
        self.ds = simulation_obj

    def run(self):
        result = [ds.current_time for ds in self.ds]
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
        result["grid_dimensions"] = self.ds.index.grid_dimensions
        result["grid_left_edges"] = self.ds.index.grid_left_edge
        result["grid_right_edges"] = self.ds.index.grid_right_edge
        result["grid_levels"] = self.ds.index.grid_levels
        result["grid_particle_count"] = self.ds.index.grid_particle_count
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
        for g in self.ds.index.grids:
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
            assert(newc == oldc)

def compare_image_lists(new_result, old_result, decimals):
    fns = []
    for i in range(2):
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        fns.append(tmpname)
    num_images = len(old_result)
    assert(num_images > 0)
    for i in range(num_images):
        mpimg.imsave(fns[0], np.loads(zlib.decompress(old_result[i])))
        mpimg.imsave(fns[1], np.loads(zlib.decompress(new_result[i])))
        results = compare_images(fns[0], fns[1], 10**(-decimals))
        if results is not None:
            if os.environ.get("JENKINS_HOME") is not None:
                tempfiles = [line.strip() for line in results.split('\n')
                             if line.endswith(".png")]
                for fn in tempfiles:
                    sys.stderr.write("\n[[ATTACHMENT|{}]]".format(fn))
                sys.stderr.write('\n')
        assert_equal(results, None, results)
        for fn in fns:
            os.remove(fn)

class VRImageComparisonTest(AnswerTestingTest):
    _type_name = "VRImageComparison"
    _attrs = ('desc',)

    def __init__(self, scene, ds, desc, decimals):
        super(VRImageComparisonTest, self).__init__(None)
        self.obj_type = ('vr',)
        self.ds = ds
        self.scene = scene
        self.desc = desc
        self.decimals = decimals

    def run(self):
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        self.scene.render()
        self.scene.save(tmpname, sigma_clip=1.0)
        image = mpimg.imread(tmpname)
        os.remove(tmpname)
        return [zlib.compress(image.dumps())]

    def compare(self, new_result, old_result):
        compare_image_lists(new_result, old_result, self.decimals)

class PlotWindowAttributeTest(AnswerTestingTest):
    _type_name = "PlotWindowAttribute"
    _attrs = ('plot_type', 'plot_field', 'plot_axis', 'attr_name', 'attr_args',
              'callback_id')
    def __init__(self, ds_fn, plot_field, plot_axis, attr_name, attr_args,
                 decimals, plot_type = 'SlicePlot', callback_id = "",
                 callback_runners = None):
        super(PlotWindowAttributeTest, self).__init__(ds_fn)
        self.plot_type = plot_type
        self.plot_field = plot_field
        self.plot_axis = plot_axis
        self.plot_kwargs = {}
        self.attr_name = attr_name
        self.attr_args = attr_args
        self.decimals = decimals
        # callback_id is so that we don't have to hash the actual callbacks
        # run, but instead we call them something
        self.callback_id = callback_id
        if callback_runners is None:
            callback_runners = []
        self.callback_runners = callback_runners

    def run(self):
        plot = self.create_plot(self.ds, self.plot_type, self.plot_field,
                                self.plot_axis, self.plot_kwargs)
        for r in self.callback_runners:
            r(self, plot)
        attr = getattr(plot, self.attr_name)
        attr(*self.attr_args[0], **self.attr_args[1])
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        plot.save(name=tmpname)
        image = mpimg.imread(tmpname)
        os.remove(tmpname)
        return [zlib.compress(image.dumps())]

    def compare(self, new_result, old_result):
        compare_image_lists(new_result, old_result, self.decimals)

class PhasePlotAttributeTest(AnswerTestingTest):
    _type_name = "PhasePlotAttribute"
    _attrs = ('plot_type', 'x_field', 'y_field', 'z_field',
              'attr_name', 'attr_args')
    def __init__(self, ds_fn, x_field, y_field, z_field,
                 attr_name, attr_args, decimals, plot_type='PhasePlot'):
        super(PhasePlotAttributeTest, self).__init__(ds_fn)
        self.data_source = self.ds.all_data()
        self.plot_type = plot_type
        self.x_field = x_field
        self.y_field = y_field
        self.z_field = z_field
        self.plot_kwargs = {}
        self.attr_name = attr_name
        self.attr_args = attr_args
        self.decimals = decimals

    def create_plot(self, data_source, x_field, y_field, z_field,
                    plot_type, plot_kwargs=None):
        # plot_type should be a string
        # plot_kwargs should be a dict
        if plot_type is None:
            raise RuntimeError('Must explicitly request a plot type')
        cls = getattr(profile_plotter, plot_type, None)
        if cls is None:
            cls = getattr(particle_plots, plot_type)
        plot = cls(*(data_source, x_field, y_field, z_field), **plot_kwargs)
        return plot

    def run(self):
        plot = self.create_plot(self.data_source, self.x_field, self.y_field,
                                self.z_field, self.plot_type, self.plot_kwargs)
        attr = getattr(plot, self.attr_name)
        attr(*self.attr_args[0], **self.attr_args[1])
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        plot.save(name=tmpname)
        image = mpimg.imread(tmpname)
        os.remove(tmpname)
        return [zlib.compress(image.dumps())]

    def compare(self, new_result, old_result):
        compare_image_lists(new_result, old_result, self.decimals)

class GenericArrayTest(AnswerTestingTest):
    _type_name = "GenericArray"
    _attrs = ('array_func_name','args','kwargs')
    def __init__(self, ds_fn, array_func, args=None, kwargs=None, decimals=None):
        super(GenericArrayTest, self).__init__(ds_fn)
        self.array_func = array_func
        self.array_func_name = array_func.__name__
        self.args = args
        self.kwargs = kwargs
        self.decimals = decimals
    def run(self):
        if self.args is None:
            args = []
        else:
            args = self.args
        if self.kwargs is None:
            kwargs = {}
        else:
            kwargs = self.kwargs
        return self.array_func(*args, **kwargs)
    def compare(self, new_result, old_result):
        if not isinstance(new_result, dict):
            new_result = {'answer': new_result}
            old_result = {'answer': old_result}

        assert_equal(len(new_result), len(old_result),
                                          err_msg="Number of outputs not equal.",
                                          verbose=True)
        for k in new_result:
            if self.decimals is None:
                assert_almost_equal(new_result[k], old_result[k])
            else:
                assert_allclose_units(new_result[k], old_result[k],
                                      10**(-self.decimals))

class GenericImageTest(AnswerTestingTest):
    _type_name = "GenericImage"
    _attrs = ('image_func_name','args','kwargs')
    def __init__(self, ds_fn, image_func, decimals, args=None, kwargs=None):
        super(GenericImageTest, self).__init__(ds_fn)
        self.image_func = image_func
        self.image_func_name = image_func.__name__
        self.args = args
        self.kwargs = kwargs
        self.decimals = decimals
    def run(self):
        if self.args is None:
            args = []
        else:
            args = self.args
        if self.kwargs is None:
            kwargs = {}
        else:
            kwargs = self.kwargs
        comp_imgs = []
        tmpdir = tempfile.mkdtemp()
        image_prefix = os.path.join(tmpdir,"test_img")
        self.image_func(image_prefix, *args, **kwargs)
        imgs = sorted(glob.glob(image_prefix+"*"))
        assert(len(imgs) > 0)
        for img in imgs:
            img_data = mpimg.imread(img)
            os.remove(img)
            comp_imgs.append(zlib.compress(img_data.dumps()))
        return comp_imgs
    def compare(self, new_result, old_result):
        compare_image_lists(new_result, old_result, self.decimals)

class AxialPixelizationTest(AnswerTestingTest):
    # This test is typically used once per geometry or coordinates type.
    # Feed it a dataset, and it checks that the results of basic pixelization
    # don't change.
    _type_name = "AxialPixelization"
    _attrs = ('geometry',)
    def __init__(self, ds_fn, decimals=None):
        super(AxialPixelizationTest, self).__init__(ds_fn)
        self.decimals = decimals
        self.geometry = self.ds.coordinates.name

    def run(self):
        rv = {}
        ds = self.ds
        for i, axis in enumerate(ds.coordinates.axis_order):
            (bounds, center, display_center) = \
                    pw.get_window_parameters(axis, ds.domain_center, None, ds)
            slc = ds.slice(axis, center[i])
            xax = ds.coordinates.axis_name[ds.coordinates.x_axis[axis]]
            yax = ds.coordinates.axis_name[ds.coordinates.y_axis[axis]]
            pix_x = ds.coordinates.pixelize(axis, slc, xax, bounds, (512, 512))
            pix_y = ds.coordinates.pixelize(axis, slc, yax, bounds, (512, 512))
            # Wipe out all NaNs
            pix_x[np.isnan(pix_x)] = 0.0
            pix_y[np.isnan(pix_y)] = 0.0
            rv['%s_x' % axis] = pix_x
            rv['%s_y' % axis] = pix_y
        return rv

    def compare(self, new_result, old_result):
        assert_equal(len(new_result), len(old_result),
                                          err_msg="Number of outputs not equal.",
                                          verbose=True)
        for k in new_result:
            if self.decimals is None:
                assert_almost_equal(new_result[k], old_result[k])
            else:
                assert_allclose_units(new_result[k], old_result[k],
                                      10**(-self.decimals))


def requires_sim(sim_fn, sim_type, big_data = False, file_check = False):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data is False and big_data is True:
        return ffalse
    elif not can_run_sim(sim_fn, sim_type, file_check):
        return ffalse
    else:
        return ftrue

def requires_answer_testing():
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if AnswerTestingTest.result_storage is not None:
        return ftrue
    else:
        return ffalse

def requires_ds(ds_fn, big_data = False, file_check = False):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data is False and big_data is True:
        return ffalse
    elif not can_run_ds(ds_fn, file_check):
        return ffalse
    else:
        return ftrue

def small_patch_amr(ds_fn, fields, input_center="max", input_weight="density"):
    if not can_run_ds(ds_fn): return
    dso = [ None, ("sphere", (input_center, (0.1, 'unitary')))]
    yield GridHierarchyTest(ds_fn)
    yield ParentageRelationshipsTest(ds_fn)
    for field in fields:
        yield GridValuesTest(ds_fn, field)
        for axis in [0, 1, 2]:
            for dobj_name in dso:
                for weight_field in [None, input_weight]:
                    yield ProjectionValuesTest(
                        ds_fn, axis, field, weight_field,
                        dobj_name)
                yield FieldValuesTest(
                        ds_fn, field, dobj_name)

def big_patch_amr(ds_fn, fields, input_center="max", input_weight="density"):
    if not can_run_ds(ds_fn):
        return
    dso = [ None, ("sphere", (input_center, (0.1, 'unitary')))]
    yield GridHierarchyTest(ds_fn)
    yield ParentageRelationshipsTest(ds_fn)
    for field in fields:
        yield GridValuesTest(ds_fn, field)
        for axis in [0, 1, 2]:
            for dobj_name in dso:
                for weight_field in [None, input_weight]:
                    yield PixelizedProjectionValuesTest(
                        ds_fn, axis, field, weight_field,
                        dobj_name)


def sph_answer(ds, ds_str_repr, ds_nparticles, fields):
    if not can_run_ds(ds):
        return
    assert_equal(str(ds), ds_str_repr)
    dso = [None, ("sphere", ("c", (0.1, 'unitary')))]
    dd = ds.all_data()
    assert_equal(dd["particle_position"].shape, (ds_nparticles, 3))
    tot = sum(dd[ptype, "particle_position"].shape[0]
              for ptype in ds.particle_types if ptype != "all")
    assert_equal(tot, ds_nparticles)
    for dobj_name in dso:
        dobj = create_obj(ds, dobj_name)
        s1 = dobj["ones"].sum()
        s2 = sum(mask.sum() for block, mask in dobj.blocks)
        assert_equal(s1, s2)
        for field, weight_field in fields.items():
            if field[0] in ds.particle_types:
                particle_type = True
            else:
                particle_type = False
            for axis in [0, 1, 2]:
                if particle_type is False:
                    yield PixelizedProjectionValuesTest(
                        ds, axis, field, weight_field,
                        dobj_name)
            yield FieldValuesTest(ds, field, dobj_name,
                                  particle_type=particle_type)

def create_obj(ds, obj_type):
    # obj_type should be tuple of
    #  ( obj_name, ( args ) )
    if obj_type is None:
        return ds.all_data()
    cls = getattr(ds, obj_type[0])
    obj = cls(*obj_type[1])
    return obj

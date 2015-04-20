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
import os
import hashlib
import contextlib
import sys
from yt.extern.six.moves import cPickle, urllib
import shelve
import zlib
import tempfile
import glob

from matplotlib.testing.compare import compare_images
from nose.plugins import Plugin
from yt.testing import *
from yt.convenience import load, simulation
from yt.config import ytcfg
from yt.data_objects.static_output import Dataset
from yt.data_objects.time_series import SimulationTimeSeries
from yt.utilities.logger import disable_stream_logging
from yt.utilities.command_line import get_yt_version

import matplotlib.image as mpimg
import yt.visualization.plot_window as pw
import yt.extern.progressbar as progressbar

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
            # Fix up filename for local storage
            if self.compare_name is not None:
                self.compare_name = "%s/%s/%s" % \
                    (os.path.realpath(options.output_dir), self.compare_name,
                     self.compare_name)
            if self.store_name is not None:
                name_dir_path = "%s/%s" % \
                    (os.path.realpath(options.output_dir),
                    self.store_name)
                if not os.path.isdir(name_dir_path):
                    os.makedirs(name_dir_path)
                self.store_name= "%s/%s" % \
                        (name_dir_path, self.store_name)
        else:
            storage_class = AnswerTestCloudStorage

        # Initialize answer/reference storage
        AnswerTestingTest.reference_storage = self.storage = \
                storage_class(self.compare_name, self.store_name)

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
        except urllib.error.HTTPError as ex:
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
        pyrax.set_credential_file(os.path.expanduser("~/.yt/rackspace"))
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
        if self.answer_name is None: return
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
    if isinstance(ds_fn, Dataset):
        return AnswerTestingTest.result_storage is not None
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    with temp_cwd(path):
        if file_check:
            return os.path.isfile(ds_fn) and \
                AnswerTestingTest.result_storage is not None
        try:
            load(ds_fn)
        except YTOutputNotIdentified:
            return False
    return AnswerTestingTest.result_storage is not None

def can_run_sim(sim_fn, sim_type, file_check = False):
    if isinstance(sim_fn, SimulationTimeSeries):
        return AnswerTestingTest.result_storage is not None
    path = ytcfg.get("yt", "test_data_dir")
    if not os.path.isdir(path):
        return False
    with temp_cwd(path):
        if file_check:
            return os.path.isfile(sim_fn) and \
                AnswerTestingTest.result_storage is not None
        try:
            simulation(sim_fn, sim_type)
        except YTOutputNotIdentified:
            return False
    return AnswerTestingTest.result_storage is not None

def data_dir_load(ds_fn, cls = None, args = None, kwargs = None):
    path = ytcfg.get("yt", "test_data_dir")
    if isinstance(ds_fn, Dataset): return ds_fn
    if not os.path.isdir(path):
        return False
    with temp_cwd(path):
        if cls is None:
            ds = load(ds_fn)
        else:
            args = args or ()
            kwargs = kwargs or {}
            ds = cls(ds_fn, *args, **kwargs)
        ds.index
        return ds

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
    def __init__(self, ds_fn):
        if isinstance(ds_fn, Dataset):
            self.ds = ds_fn
        else:
            self.ds = data_dir_load(ds_fn)

    def __call__(self):
        nv = self.run()
        if self.reference_storage.reference_name is not None:
            dd = self.reference_storage.get(self.storage_name)
            if dd is None or self.description not in dd:
                raise YTNoOldAnswer("%s : %s" % (self.storage_name , self.description))
            ov = dd[self.description]
            self.compare(nv, ov)
        else:
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
        # plot_args should be a tuple
        # plot_kwargs should be a dict
        if plot_type is None:
            raise RuntimeError('Must explicitly request a plot type')
        cls = getattr(pw, plot_type)
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
        return "_".join(args)

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
        if self.particle_type:
            weight_field = "particle_ones"
        else:
            weight_field = "ones"
        avg = obj.quantities.weighted_average_quantity(
            self.field, weight=weight_field)
        mi, ma = obj.quantities.extrema(self.field)
        return np.array([avg, mi, ma])

    def compare(self, new_result, old_result):
        err_msg = "Field values for %s not equal." % (self.field,)
        if self.decimals is None:
            assert_equal(new_result, old_result,
                         err_msg=err_msg, verbose=True)
        else:
            assert_allclose(new_result, old_result, 10.**(-self.decimals),
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
                assert_allclose(nres, ores, 10.**-(self.decimals),
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

class SimulatedHaloMassFunctionTest(AnswerTestingTest):
    _type_name = "SimulatedHaloMassFunction"
    _attrs = ("finder",)

    def __init__(self, ds_fn, finder):
        super(SimulatedHaloMassFunctionTest, self).__init__(ds_fn)
        self.finder = finder
    
    def run(self):
        from yt.analysis_modules.halo_analysis.api import HaloCatalog
        from yt.analysis_modules.halo_mass_function.api import HaloMassFcn
        hc = HaloCatalog(data_ds=self.ds, finder_method=self.finder)
        hc.create()
        
        hmf = HaloMassFcn(halos_ds=hc.halos_ds)
        result = np.empty((2, hmf.masses_sim.size))
        result[0] = hmf.masses_sim.d
        result[1] = hmf.n_cumulative_sim.d
        return result

    def compare(self, new_result, old_result):
        err_msg = ("Simulated halo mass functions not equation for " +
                   "%s halo finder.") % self.finder
        assert_equal(new_result, old_result,
                     err_msg=err_msg, verbose=True)

class AnalyticHaloMassFunctionTest(AnswerTestingTest):
    _type_name = "AnalyticHaloMassFunction"
    _attrs = ("fitting_function",)

    def __init__(self, ds_fn, fitting_function):
        super(AnalyticHaloMassFunctionTest, self).__init__(ds_fn)
        self.fitting_function = fitting_function
    
    def run(self):
        from yt.analysis_modules.halo_mass_function.api import HaloMassFcn
        hmf = HaloMassFcn(simulation_ds=self.ds,
                          fitting_function=self.fitting_function)
        result = np.empty((2, hmf.masses_analytic.size))
        result[0] = hmf.masses_analytic.d
        result[1] = hmf.n_cumulative_analytic.d
        return result

    def compare(self, new_result, old_result):
        err_msg = ("Analytic halo mass functions not equation for " +
                   "fitting function %d.") % self.fitting_function
        assert_equal(new_result, old_result,
                     err_msg=err_msg, verbose=True)

def compare_image_lists(new_result, old_result, decimals):
    fns = ['old.png', 'new.png']
    num_images = len(old_result)
    assert(num_images > 0)
    for i in range(num_images):
        mpimg.imsave(fns[0], np.loads(zlib.decompress(old_result[i])))
        mpimg.imsave(fns[1], np.loads(zlib.decompress(new_result[i])))
        assert compare_images(fns[0], fns[1], 10**(-decimals)) == None
        for fn in fns: os.remove(fn)

class PlotWindowAttributeTest(AnswerTestingTest):
    _type_name = "PlotWindowAttribute"
    _attrs = ('plot_type', 'plot_field', 'plot_axis', 'attr_name', 'attr_args')
    def __init__(self, ds_fn, plot_field, plot_axis, attr_name, attr_args,
                 decimals, plot_type = 'SlicePlot'):
        super(PlotWindowAttributeTest, self).__init__(ds_fn)
        self.plot_type = plot_type
        self.plot_field = plot_field
        self.plot_axis = plot_axis
        self.plot_kwargs = {}
        self.attr_name = attr_name
        self.attr_args = attr_args
        self.decimals = decimals

    def run(self):
        plot = self.create_plot(self.ds, self.plot_type, self.plot_field,
                                self.plot_axis, self.plot_kwargs)
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
        assert_equal(len(new_result), len(old_result),
                                          err_msg="Number of outputs not equal.",
                                          verbose=True)
        for k in new_result:
            if self.decimals is None:
                assert_equal(new_result[k], old_result[k])
            else:
                assert_allclose(new_result[k], old_result[k], 10**(-self.decimals))

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
        imgs = glob.glob(image_prefix+"*")
        assert(len(imgs) > 0)
        for img in imgs:
            img_data = mpimg.imread(img)
            os.remove(img)
            comp_imgs.append(zlib.compress(img_data.dumps()))
        return comp_imgs
    def compare(self, new_result, old_result):
        compare_image_lists(new_result, old_result, self.decimals)

def requires_sim(sim_fn, sim_type, big_data = False, file_check = False):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data == False and big_data == True:
        return ffalse
    elif not can_run_sim(sim_fn, sim_type, file_check):
        return ffalse
    else:
        return ftrue
        
def requires_ds(ds_fn, big_data = False, file_check = False):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if run_big_data == False and big_data == True:
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
    if not can_run_ds(ds_fn): return
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

def create_obj(ds, obj_type):
    # obj_type should be tuple of
    #  ( obj_name, ( args ) )
    if obj_type is None:
        return ds.all_data()
    cls = getattr(ds, obj_type[0])
    obj = cls(*obj_type[1])
    return obj

class AssertWrapper(object):
    """
    Used to wrap a numpy testing assertion, in order to provide a useful name
    for a given assertion test.
    """
    def __init__(self, description, *args):
        # The key here is to add a description attribute, which nose will pick
        # up.
        self.args = args
        self.description = description

    def __call__(self):
        self.args[0](*self.args[1:])

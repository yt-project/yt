"""
Answer Testing using pytest as a starting point



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
from yt.funcs import \
    get_pbar
from yt.testing import \
    assert_allclose_units
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

@contextlib.contextmanager
def temp_cwd(cwd):
    oldcwd = os.getcwd()
    os.chdir(cwd)
    yield
    os.chdir(oldcwd)

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

class FieldValuesTest():
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

class AllFieldValuesTest():
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

class ProjectionValuesTest():
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

class PixelizedProjectionValuesTest():
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

class GridValuesTest():
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

class VerifySimulationSameTest():
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

class GridHierarchyTest():
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

class ParentageRelationshipsTest():
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

class SimulatedHaloMassFunctionTest():
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

class AnalyticHaloMassFunctionTest():
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
        err_msg = ("Analytic halo mass functions not equal for " +
                   "fitting function %d.") % self.fitting_function
        assert_almost_equal(new_result, old_result,
                            err_msg=err_msg, verbose=True)

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

class VRImageComparisonTest():
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

class PlotWindowAttributeTest():
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

class PhasePlotAttributeTest():
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

class GenericArrayTest():
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

class GenericImageTest():
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

class AxialPixelizationTest():
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

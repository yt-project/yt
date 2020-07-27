"""
Title: answer_tests.py
Purpose: Contains answer tests that are used by yt's various frontends
"""
import glob
import os
import tempfile

import matplotlib.image as mpimg
import numpy as np

from yt.testing import assert_equal
from yt.utilities.on_demand_imports import _h5py as h5py
import yt.visualization.particle_plots as particle_plots
import yt.visualization.plot_window as pw

from . import utils


def grid_hierarchy(ds):
    result = {}
    result["grid_dimensions"] = ds.index.grid_dimensions
    result["grid_left_edges"] = ds.index.grid_left_edge
    result["grid_right_edges"] = ds.index.grid_right_edge
    result["grid_levels"] = ds.index.grid_levels
    result["grid_particle_count"] = ds.index.grid_particle_count
    return result


def parentage_relationships(ds):
    result = {}
    result["parents"] = []
    result["children"] = []
    for g in ds.index.grids:
        p = g.Parent
        if p is None:
            result["parents"].append(-1)
        elif hasattr(p, "id"):
            result["parents"].append(p.id)
        else:
            result["parents"].append([pg.id for pg in p])
        result["children"].append([c.id for c in g.Children])
    result["parents"] = np.array(result["parents"])
    result["children"] = np.array(result["children"])
    return result


def grid_values(ds, field):
    result = {}
    for g in ds.index.grids:
        result[str(g.id)] = g[field]
        g.clear_data()
    return result


def projection_values(ds, axis, field, weight_field=None, dobj_type=None):
    if dobj_type is not None:
        dobj = utils.create_obj(ds, dobj_type)
    else:
        dobj = None
    if ds.domain_dimensions[axis] == 1:
        return None
    proj = ds.proj(field, axis, weight_field=weight_field, data_source=dobj)
    return proj.field_data


def field_values(ds, field, obj_type=None, particle_type=False):
    # If needed build an instance of the dataset type
    obj = utils.create_obj(ds, obj_type)
    determined_field = obj._determine_fields(field)[0]
    fd = ds.field_info[field]
    # Get the proper weight field depending on if we're looking at
    # particles or not
    if particle_type:
        weight_field = (determined_field[0], "particle_ones")
    elif fd.is_sph_field:
        weight_field = (determined_field[0], "ones")
    else:
        weight_field = ("index", "ones")
    # Get the average, min, and max
    avg = obj.quantities.weighted_average_quantity(
        determined_field, weight=weight_field
    )
    minimum, maximum = obj.quantities.extrema(field)
    # Return as a hashable bytestring
    return np.array([avg, minimum, maximum])


def all_field_values(ds, field, obj_type=None):
    obj = utils.create_obj(ds, obj_type)
    return obj[field]


def pixelized_projection_values(ds, axis, field, weight_field=None, dobj_type=None):
    if dobj_type is not None:
        obj = utils.create_obj(ds, dobj_type)
    else:
        obj = None
    proj = ds.proj(field, axis, weight_field=weight_field, data_source=obj)
    frb = proj.to_frb((1.0, "unitary"), 256)
    frb[field]
    if weight_field is not None:
        frb[weight_field]
    d = frb.data
    for f in proj.field_data:
        # Sometimes f will be a tuple.
        d["%s_sum" % (f,)] = proj.field_data[f].sum(dtype="float64")
    return d


def pixelized_particle_projection_values(
    ds, axis, field, weight_field=None, dobj_type=None
):
    proj_plot = particle_plots.ParticleProjectionPlot(
        ds, axis, [field], weight_field=weight_field
    )
    proj = proj_plot.data_source
    frb = proj_plot.frb
    frb[field]
    if weight_field is not None:
        frb[weight_field]
    d = frb.data
    for f in proj.field_data:
        # Sometimes f will be a tuple.
        d["%s_sum" % (f,)] = proj.field_data[f].sum(dtype="float64")
    return d


def small_patch_amr(ds, field, weight, axis, ds_obj):
    results = {}
    results["grid_hierarchy"] = grid_hierarchy(ds)
    results["parentage_relationships"] = parentage_relationships(ds)
    results["grid_values"] = grid_values(ds, field)
    results["field_values"] = field_values(ds, field, ds_obj)
    results["projection_values"] = projection_values(ds, axis, field, weight, ds_obj)
    return results


def big_patch_amr(ds, field, weight, axis, ds_obj):
    results = {}
    results["grid_hierarchy"] = grid_hierarchy(ds)
    results["parentage_relationships"] = parentage_relationships(ds)
    results["grid_values"] = grid_values(ds, field)
    results["pixelized_projection_values"] = pixelized_projection_values(
        ds, axis, field, weight, ds_obj
    )
    return results


def generic_array(func, args=[], kwargs={}):
    return func(*args, **kwargs)


def sph_answer(ds, ds_str_repr, ds_nparticles, field, weight, ds_obj, axis):
    assert str(ds) == ds_str_repr
    hex_digests = {}
    dd = ds.all_data()
    assert_equal(dd["all", "particle_position"].shape, (ds_nparticles, 3))
    tot = sum(
        dd[ptype, "particle_position"].shape[0] for ptype in ds.particle_types_raw
    )
    assert_equal(tot, ds_nparticles)
    particle_type = field[0] in ds.particle_types
    if not particle_type:
        ppv = pixelized_projection_values(ds, axis, field, weight, ds_obj)
        hex_digests["pixelized_projection_values"] = ppv
    fv = field_values(ds, field, ds_obj, particle_type=particle_type)
    hex_digests["field_values"] = fv
    return hex_digests


def nbody_answer(ds, ds_str_repr, ds_nparticles, field, weight, ds_obj, axis):
    assert str(ds) == ds_str_repr
    hex_digests = {}
    dd = ds.all_data()
    assert_equal(dd["all", "particle_position"].shape, (ds_nparticles, 3))
    tot = sum(
        dd[ptype, "particle_position"].shape[0] for ptype in ds.particle_types_raw
    )
    assert_equal(tot, ds_nparticles)
    particle_type = field[0] in ds.particle_types
    if not particle_type:
        pppv = pixelized_particle_projection_values(ds, axis, field, weight, ds_obj)
        hex_digests["pixelized_particle_projection_values"] = pppv
    fv = field_values(ds, field, ds_obj, particle_type=particle_type)
    hex_digests["field_values"] = fv
    return hex_digests


def get_field_size_and_mean(ds, field, geometric):
    if geometric:
        obj = ds.all_data()
    else:
        obj = ds.data
    return np.array([obj[field].size, obj[field].mean()])


def plot_window_attribute(
    ds,
    plot_field,
    plot_axis,
    attr_name,
    attr_args,
    plot_type="SlicePlot",
    callback_id="",
    callback_runners=None,
):
    if callback_runners is None:
        callback_runners = []
    plot = utils._create_plot_window_attribute_plot(
        ds, plot_type, plot_field, plot_axis, {}
    )
    for r in callback_runners:
        r(plot_field, plot)
    attr = getattr(plot, attr_name)
    attr(*attr_args[0], **attr_args[1])
    tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
    os.close(tmpfd)
    plot.save(name=tmpname)
    image = mpimg.imread(tmpname)
    os.remove(tmpname)
    return image


def phase_plot_attribute(
    ds,
    x_field,
    y_field,
    z_field,
    attr_name,
    attr_args,
    plot_type="PhasePlot",
    plot_kwargs={},
):
    data_source = ds.all_data()
    plot = utils._create_phase_plot_attribute_plot(
        data_source, x_field, y_field, z_field, plot_type, plot_kwargs
    )
    attr = getattr(plot, attr_name)
    attr(*attr_args[0], **attr_args[1])
    tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
    os.close(tmpfd)
    plot.save(name=tmpname)
    image = mpimg.imread(tmpname)
    os.remove(tmpname)
    return image


def generic_image(img_func, args=[], kwargs={}):
    comp_imgs = []
    tmpdir = tempfile.mkdtemp()
    image_prefix = os.path.join(tmpdir, "test_img")
    img_func(image_prefix, *args, **kwargs)
    imgs = sorted(glob.glob(image_prefix + "*"))
    assert len(imgs) > 0
    for img in imgs:
        img_data = mpimg.imread(img)
        os.remove(img)
        comp_imgs.append(img_data)
    return comp_imgs


def axial_pixelization(ds):
    r"""
    This test is typically used once per geometry or coordinates type.
    Feed it a dataset, and it checks that the results of basic pixelization
    don't change.
    """
    rv = {}
    for i, axis in enumerate(ds.coordinates.axis_order):
        (bounds, center, display_center) = pw.get_window_parameters(
            axis, ds.domain_center, None, ds
        )
        slc = ds.slice(axis, center[i])
        xax = ds.coordinates.axis_name[ds.coordinates.x_axis[axis]]
        yax = ds.coordinates.axis_name[ds.coordinates.y_axis[axis]]
        pix_x = ds.coordinates.pixelize(axis, slc, xax, bounds, (512, 512))
        pix_y = ds.coordinates.pixelize(axis, slc, yax, bounds, (512, 512))
        # Wipe out invalid values (fillers)
        pix_x[~np.isfinite(pix_x)] = 0.0
        pix_y[~np.isfinite(pix_y)] = 0.0
        rv["%s_x" % axis] = pix_x
        rv["%s_y" % axis] = pix_y
    return rv


def extract_connected_sets(ds_fn, data_source, field, num_levels, min_val, max_val):
    n, all_sets = data_source.extract_connected_sets(
        field, num_levels, min_val, max_val
    )
    result = []
    for level in all_sets:
        for set_id in all_sets[level]:
            result.append(
                [
                    all_sets[level][set_id]["cell_mass"].size,
                    all_sets[level][set_id]["cell_mass"].sum(),
                ]
            )
    result = np.array(result)
    return result


def VR_image_comparison(scene):
    tmpfd, tmpname = tempfile.mkstemp(suffix=".png")
    os.close(tmpfd)
    scene.save(tmpname, sigma_clip=1.0)
    image = mpimg.imread(tmpname)
    os.remove(tmpname)
    return image


def yt_data_field(ds, field, geometric):
    if geometric:
        obj = ds.all_data()
    else:
        obj = ds.data
    num_e = obj[field].size
    avg = obj[field].mean()
    return np.array([num_e, avg])


def verify_simulation_same(sim_obj):
    result = [ds.current_time for ds in sim_obj]
    return np.array(result)

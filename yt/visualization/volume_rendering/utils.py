import numpy as np

from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer3D,
)
from yt.data_objects.static_output import Dataset
from yt.utilities.lib import bounding_volume_hierarchy
from yt.utilities.lib.image_samplers import (
    InterpolatedProjectionSampler,
    ProjectionSampler,
    VolumeRenderSampler,
)
from yt.utilities.on_demand_imports import NotAModule

try:
    from yt.utilities.lib.embree_mesh import mesh_traversal  # type: ignore
# Catch ValueError in case size of objects in Cython change
except (ImportError, ValueError):
    mesh_traversal = NotAModule("pyembree")


def data_source_or_all(data_source):
    if isinstance(data_source, Dataset):
        data_source = data_source.all_data()
    if not isinstance(data_source, (YTSelectionContainer3D, type(None))):
        raise RuntimeError(
            "The data_source is not a valid 3D data container.\n"
            "Expected an object of type YTSelectionContainer3D but received "
            "an object of type %s." % type(data_source)
        )
    return data_source


def new_mesh_sampler(camera, render_source, engine):
    params = ensure_code_unit_params(camera._get_sampler_params(render_source))
    args = (
        np.atleast_3d(params["vp_pos"]),
        np.atleast_3d(params["vp_dir"]),
        params["center"],
        params["bounds"],
        np.atleast_3d(params["image"]).astype("float64"),
        params["x_vec"],
        params["y_vec"],
        params["width"],
        render_source.volume_method,
    )
    kwargs = {"lens_type": params["lens_type"]}
    if engine == "embree":
        sampler = mesh_traversal.EmbreeMeshSampler(*args, **kwargs)
    elif engine == "yt":
        sampler = bounding_volume_hierarchy.BVHMeshSampler(*args, **kwargs)
    return sampler


def new_volume_render_sampler(camera, render_source):
    params = ensure_code_unit_params(camera._get_sampler_params(render_source))
    params.update(transfer_function=render_source.transfer_function)
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        np.atleast_3d(params["vp_pos"]),
        np.atleast_3d(params["vp_dir"]),
        params["center"],
        params["bounds"],
        params["image"],
        params["x_vec"],
        params["y_vec"],
        params["width"],
        render_source.volume_method,
        params["transfer_function"],
        params["num_samples"],
    )
    kwargs = {
        "lens_type": params["lens_type"],
    }
    if "camera_data" in params:
        kwargs["camera_data"] = params["camera_data"]
    if render_source.zbuffer is not None:
        kwargs["zbuffer"] = render_source.zbuffer.z
        args[4][:] = np.reshape(
            render_source.zbuffer.rgba[:],
            (camera.resolution[0], camera.resolution[1], 4),
        )
    else:
        kwargs["zbuffer"] = np.ones(params["image"].shape[:2], "float64")
    sampler = VolumeRenderSampler(*args, **kwargs)
    return sampler


def new_interpolated_projection_sampler(camera, render_source):
    params = ensure_code_unit_params(camera._get_sampler_params(render_source))
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        np.atleast_3d(params["vp_pos"]),
        np.atleast_3d(params["vp_dir"]),
        params["center"],
        params["bounds"],
        params["image"],
        params["x_vec"],
        params["y_vec"],
        params["width"],
        render_source.volume_method,
        params["num_samples"],
    )
    kwargs = {"lens_type": params["lens_type"]}
    if render_source.zbuffer is not None:
        kwargs["zbuffer"] = render_source.zbuffer.z
    else:
        kwargs["zbuffer"] = np.ones(params["image"].shape[:2], "float64")
    sampler = InterpolatedProjectionSampler(*args, **kwargs)
    return sampler


def new_projection_sampler(camera, render_source):
    params = ensure_code_unit_params(camera._get_sampler_params(render_source))
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        np.atleast_3d(params["vp_pos"]),
        np.atleast_3d(params["vp_dir"]),
        params["center"],
        params["bounds"],
        params["image"],
        params["x_vec"],
        params["y_vec"],
        params["width"],
        render_source.volume_method,
        params["num_samples"],
    )
    kwargs = {
        "lens_type": params["lens_type"],
    }
    if render_source.zbuffer is not None:
        kwargs["zbuffer"] = render_source.zbuffer.z
    else:
        kwargs["zbuffer"] = np.ones(params["image"].shape[:2], "float64")
    sampler = ProjectionSampler(*args, **kwargs)
    return sampler


def get_corners(le, re):
    return np.array(
        [
            [le[0], le[1], le[2]],
            [re[0], le[1], le[2]],
            [re[0], re[1], le[2]],
            [le[0], re[1], le[2]],
            [le[0], le[1], re[2]],
            [re[0], le[1], re[2]],
            [re[0], re[1], re[2]],
            [le[0], re[1], re[2]],
        ],
        dtype="float64",
    )


def ensure_code_unit_params(params):
    for param_name in ["center", "vp_pos", "vp_dir", "width"]:
        param = params[param_name]
        if hasattr(param, "in_units"):
            params[param_name] = param.in_units("code_length")
    bounds = params["bounds"]
    if hasattr(bounds[0], "units"):
        params["bounds"] = tuple(b.in_units("code_length").d for b in bounds)
    return params

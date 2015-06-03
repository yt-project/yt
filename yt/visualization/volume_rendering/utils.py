import numpy as np
from yt.data_objects.static_output import Dataset
from yt.utilities.lib.grid_traversal import \
    VolumeRenderSampler, InterpolatedProjectionSampler, ProjectionSampler
from yt.utilities.lib.mesh_traversal import \
    MeshSampler


def data_source_or_all(data_source):
    if isinstance(data_source, Dataset):
        data_source = data_source.all_data()
    return data_source


def new_mesh_render_sampler(camera, render_source):
    params = camera._get_sampler_params(render_source)
    args = (
        params['vp_pos'],
        params['vp_dir'],
        params['center'],
        params['bounds'],
        params['image'],
        params['x_vec'],
        params['y_vec'],
        params['width'],
    )

    sampler = MeshSampler(*args)
    return sampler


def new_volume_render_sampler(camera, render_source):
    params = camera._get_sampler_params(render_source)
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        params['vp_pos'],
        params['vp_dir'],
        params['center'],
        params['bounds'],
        params['image'],
        params['x_vec'],
        params['y_vec'],
        params['width'],
        params['transfer_function'],
        params['num_samples'],
    )
    kwargs = {}
    if render_source.zbuffer is not None:
        kwargs['zbuffer'] = render_source.zbuffer.z
        args[4][:] = render_source.zbuffer.rgba[:]

    sampler = VolumeRenderSampler(*args, **kwargs)
    return sampler


def new_interpolated_projection_sampler(camera, render_source):
    params = camera._get_sampler_params(render_source)
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        params['vp_pos'],
        params['vp_dir'],
        params['center'],
        params['bounds'],
        params['image'],
        params['x_vec'],
        params['y_vec'],
        params['width'],
        params['num_samples'],
    )
    kwargs = {}
    if render_source.zbuffer is not None:
        kwargs['zbuffer'] = render_source.zbuffer.z
    sampler = InterpolatedProjectionSampler(*args)
    return sampler


def new_projection_sampler(camera, render_source):
    params = camera._get_sampler_params(render_source)
    params.update(transfer_function=render_source.transfer_function)
    params.update(num_samples=render_source.num_samples)
    args = (
        params['vp_pos'],
        params['vp_dir'],
        params['center'],
        params['bounds'],
        params['image'],
        params['x_vec'],
        params['y_vec'],
        params['width'],
        params['num_samples'],
    )
    kwargs = {}
    if render_source.zbuffer is not None:
        kwargs['zbuffer'] = render_source.zbuffer.z
    sampler = ProjectionSampler(*args)
    return sampler


def get_corners(le, re):
    return np.array([
        [le[0], le[1], le[2]],
        [re[0], le[1], le[2]],
        [re[0], re[1], le[2]],
        [le[0], re[1], le[2]],
        [le[0], le[1], re[2]],
        [re[0], le[1], re[2]],
        [re[0], re[1], re[2]],
        [le[0], re[1], re[2]],
        ], dtype='float64')


from yt.data_objects.static_output import Dataset
from yt.utilities.lib.grid_traversal import \
    VolumeRenderSampler


def data_source_or_all(data_source):
    if isinstance(data_source, Dataset):
        data_source = data_source.all_data()
    return data_source


def new_volume_render_sampler(camera, render_source):
    params = camera.get_sampler_params(render_source)
    params.update(transfer_function=render_source.transfer_function)
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
    sampler = VolumeRenderSampler(*args, **kwargs)
    return sampler

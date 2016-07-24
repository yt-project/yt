import numpy as np

from yt.testing import \
    fake_random_ds, \
    assert_array_equal
from yt.units.yt_array import \
    YTArray, \
    uintersect1d

def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"

# Copied from test_boolean for computing a unique identifier for
# each cell from cell positions
def _IDFIELD(field, data):
    width = data.ds.domain_right_edge - data.ds.domain_left_edge
    min_dx = YTArray(1.0/8192, input_units='code_length',
                     registry=data.ds.unit_registry)
    delta = width / min_dx
    x = data['x'] - min_dx / 2.
    y = data['y'] - min_dx / 2.
    z = data['z'] - min_dx / 2.
    xi = x / min_dx
    yi = y / min_dx
    zi = z / min_dx
    index = xi + delta[0] * (yi + delta[1] * zi)
    return index

def test_compose_no_overlap():
    r"""Test to make sure that composed data objects that don't
    overlap behave the way we expect (return empty collections)
    """
    empty = np.array([])
    for n in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs=n)
        ds.add_field(("index", "ID"), function=_IDFIELD)

        # position parameters for initial region
        center = [0.25]*3
        left_edge = [0.1]*3
        right_edge = [0.4]*3
        normal = [1, 0, 0]
        radius = height = 0.15

        # initial 3D regions
        sources = [ds.sphere(center, radius),
                   ds.region(center, left_edge, right_edge),
                   ds.disk(center, normal, radius, height)]

        # position parameters for non-overlapping regions
        center = [0.75]*3
        left_edge = [0.6]*3
        right_edge = [0.9]*3

        # subselect non-overlapping 0, 1, 2, 3D regions
        for data1 in sources:
            data2 = ds.sphere(center, radius, data_source=data1)
            yield assert_array_equal, data2['index', 'ID'], empty

            data2 = ds.region(center, left_edge, right_edge, data_source=data1)
            yield assert_array_equal, data2['index', 'ID'], empty  

            data2 = ds.disk(center, normal, radius, height, data_source=data1)
            yield assert_array_equal, data2['index', 'ID'], empty

            for d in range(3):
                data2 = ds.slice(d, center[d], data_source=data1)
                yield assert_array_equal, data2['index', 'ID'], empty

            for d in range(3):
                data2 = ds.ortho_ray(d, center[0:d] + center[d+1:], data_source=data1)
                yield assert_array_equal, data2['index', 'ID'], empty

            data2 = ds.point(center, data_source=data1)
            yield assert_array_equal, data2['index', 'ID'], empty

def test_compose_overlap():
    r"""Test to make sure that composed data objects that do
    overlap behave the way we expect
    """
    for n in [1, 2, 4, 8]:
        ds = fake_random_ds(64, nprocs=n)
        ds.add_field(("index", "ID"), function=_IDFIELD)

        # position parameters for initial region
        center = [0.4, 0.5, 0.5]
        left_edge = [0.1]*3
        right_edge = [0.7]*3
        normal = [1, 0, 0]
        radius = height = 0.15

        # initial 3D regions
        sources = [ds.sphere(center, radius),
                   ds.region(center, left_edge, right_edge),
                   ds.disk(center, normal, radius, height)]

        # position parameters for overlapping regions
        center = [0.6, 0.5, 0.5]
        left_edge = [0.3]*3
        right_edge = [0.9]*3

        # subselect non-overlapping 0, 1, 2, 3D regions
        for data1 in sources:
            id1 = data1['index', 'ID']

            data2 = ds.sphere(center, radius)
            data3 = ds.sphere(center, radius, data_source=data1)
            id2 = data2['index', 'ID']
            id3 = data3['index', 'ID']
            id3.sort()
            yield assert_array_equal, uintersect1d(id1, id2), id3

            data2 = ds.region(center, left_edge, right_edge)
            data3 = ds.region(center, left_edge, right_edge, data_source=data1)
            id2 = data2['index', 'ID']
            id3 = data3['index', 'ID']
            id3.sort()
            yield assert_array_equal, uintersect1d(id1, id2), id3

            data2 = ds.disk(center, normal, radius, height)
            data3 = ds.disk(center, normal, radius, height, data_source=data1)
            id2 = data2['index', 'ID']
            id3 = data3['index', 'ID']
            id3.sort()
            yield assert_array_equal, uintersect1d(id1, id2), id3

            for d in range(3):
                data2 = ds.slice(d, center[d])
                data3 = ds.slice(d, center[d], data_source=data1)
                id2 = data2['index', 'ID']
                id3 = data3['index', 'ID']
                id3.sort()
                yield assert_array_equal, uintersect1d(id1, id2), id3

            for d in range(3):
                data2 = ds.ortho_ray(d, center[0:d] + center[d+1:])
                data3 = ds.ortho_ray(d, center[0:d] + center[d+1:], data_source=data1)
                id2 = data2['index', 'ID']
                id3 = data3['index', 'ID']
                id3.sort()
                yield assert_array_equal, uintersect1d(id1, id2), id3

            data2 = ds.point(center)
            data3 = ds.point(center, data_source=data1)
            id2 = data2['index', 'ID']
            id3 = data3['index', 'ID']
            id3.sort()
            yield assert_array_equal, uintersect1d(id1, id2), id3

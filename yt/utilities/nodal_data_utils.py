import numpy as np

_index_map = np.array(
    [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 1, 0, 1],
        [0, 0, 1, 1, 0, 0, 1, 1],
        [0, 1, 2, 3, 0, 1, 2, 3],
        [0, 0, 0, 0, 1, 1, 1, 1],
        [0, 1, 0, 1, 2, 3, 2, 3],
        [0, 0, 1, 1, 2, 2, 3, 3],
        [0, 1, 2, 3, 4, 5, 6, 7],
    ]
)


def _get_linear_index(nodal_flag):
    if len(nodal_flag) == 2:
        return 1 * nodal_flag[1] + 2 * nodal_flag[0]
    else:
        return 1 * nodal_flag[2] + 2 * nodal_flag[1] + 4 * nodal_flag[0]


def _get_indices(nodal_flag):
    li = _get_linear_index(nodal_flag)
    return _index_map[li]


def get_nodal_data(data_source, field):
    finfo = data_source.ds._get_field_info(field)
    nodal_flag = finfo.nodal_flag
    field_data = data_source[field]
    inds = _get_indices(nodal_flag)
    return field_data[:, inds]


def get_nodal_slices(shape, nodal_flag, dim):
    slices = []
    dir_slices = [[] for _ in range(dim)]

    for i in range(dim):
        if nodal_flag[i]:
            dir_slices[i] = [slice(0, shape[i] - 1), slice(1, shape[i])]
        else:
            dir_slices[i] = [slice(0, shape[i])]

    for sl_i in dir_slices[0]:
        for sl_j in dir_slices[1]:
            if dim > 2:
                for sl_k in dir_slices[2]:
                    slices.append([sl_i, sl_j, sl_k])
            else:
                slices.append([sl_i, sl_j])

    return tuple(slices)

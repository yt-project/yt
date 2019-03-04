import yt
import numpy as np

ds = yt.load('output_00080/info_00080.txt')
sp = ds.sphere([.5]*3, (0.5, 'code_length'))
ad = ds.all_data()

# def simple(field ,data):
#     if isinstance(data, yt.fields.field_detector.FieldDetector):
#         return data['pressure'] / data['dx']

#     # dest = np.zeros((6_000_000), dtype=np.float64) * np.nan
#     offset = 0
#     cell_count = 0
#     import ipdb; ipdb.set_trace()

#     for i, subset in enumerate(data._current_chunk.objs):
#         oh = subset.domain.oct_handler
#         cell_count += data.selector.count_oct_cells(oh, subset.domain_id)
    
#     dest = np.zeros(cell_count, dtype=np.float64)

#     for i, subset in enumerate(data._current_chunk.objs):
#         # Extract *all* data in octree
#         tmp = subset['x']
#         doffset = subset.select(subset.selector, tmp, dest, offset)
#         offset += doffset

#     return data.apply_units(dest, tmp.units)

def simple(field ,data):
    if isinstance(data, yt.fields.field_detector.FieldDetector):
        return data['pressure'] / data['dx']

    chunks = list(data.index._chunk_io(data))

    all_data = []
    for ichunk, chunk in enumerate(chunks):
        for subset in chunk.objs:
            tmp = subset['x'].T.reshape(-1, 8)# .T.reshape(-1, 8)
            oh = subset.domain.oct_handler
            selector = data.selector
            dom_cell_count = data.selector.count_oct_cells(oh, subset.domain_id)
            print(np.product(tmp.shape), dom_cell_count)
            levels, cell_inds, file_inds = oh.file_index_octs(
                selector, subset.domain_id, dom_cell_count)

            tr = {}
            tr['x'] = np.zeros(dom_cell_count, 'float64')
            for ilevel in range(levels.max()):
                oh.fill_level(ilevel, levels, cell_inds, file_inds, tr, {'x': tmp})

            all_data.append(tr['x'])

    dest = np.concatenate(all_data)

    return data.apply_units(dest, tmp.units)

ds.add_field(('gas', 'test'), function=simple, units='code_length')
test = ad['gas', 'test']
a = test.to('unitary').value
b = sp['dx'].to('unitary').value
assert np.allclose(a, b)

assert test.shape == ad['ones'].shape
print('Yes!')

import sys
sys.exit(0)

def generate_gradient(direction):
    idir = 'xyz'.index(direction)
    def grad(field, data):
        if isinstance(data, yt.fields.field_detector.FieldDetector):
            return data['pressure'] / data['dx']
        offset = 0
        cell_count = 0

        for i, subset in enumerate(data._current_chunk.objs):
            oh = subset.domain.oct_handler
            cell_count += data.selector.count_oct_cells(oh, subset.domain_id)
        
        dest = np.zeros(cell_count, dtype=np.float64)

        for i, subset in enumerate(data._current_chunk.objs):
            oh = subset.domain.oct_handler
            # # Extract *all* data in octree
            # data_in = {direction: subset[direction],
            #            'pressure': subset['pressure']}
            # data_out = oh.get_hypercube(subset, data_in)
            # # xin = data_in['x'][..., 100].to('code_length').value*64 + .5
            # # xout = data_out['x'][1:3, 1:3, 1:3, 100].to('code_length').value*64 + .5
            # sl = slice(1, 3)
            # sl1 = [sl]*3 + [slice(None)]
            # sl2 = [sl]*3 + [slice(None)]
            # sl3 = [sl]*3 + [slice(None)]
            # sl1[2-idir] = slice(0, 2)
            # sl2[2-idir] = slice(1, 3)
            # sl3[2-idir] = slice(2, 4)

            # # Compute gradients
            # x = data_out[direction]
            # p = data_out['pressure']
            # dpl = (p[sl3] - p[sl2]) / (x[sl3] - x[sl2])
            # dpr = (p[sl2] - p[sl1]) / (x[sl2] - x[sl1])

            # maskl = np.isfinite(dpl)
            # maskr = np.isfinite(dpr)

            # grad = np.where(maskl & maskr, (dpl + dpr) / 2,
            #     np.where(maskl, dpl, dpr))

            # Select data in region
            tmp = subset['x']
            doffset = oh.selector_fill(subset.selector, tmp, dest, offset)
            offset += doffset

        return data.apply_units(dest, p.units / x.units)
    return grad
# ds.add_field(('gas', 'p_grad_x'), function=generate_gradient('x'), sampling_type='cell', units='dyne/cm**3')
# ds.add_field(('gas', 'p_grad_y'), function=generate_gradient('y'), sampling_type='cell', units='dyne/cm**3')
# ds.add_field(('gas', 'p_grad_z'), function=generate_gradient('z'), sampling_type='cell', units='dyne/cm**3')
# test = ad['gas', 'p_grad_x']
# mask = np.isfinite(test)

# p = yt.ProjectionPlot(ds, 'x', 'p_grad_x')


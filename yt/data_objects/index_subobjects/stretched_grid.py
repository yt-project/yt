import numpy as np

from yt.geometry.selection_routines import convert_mask_to_indices

from .grid_patch import AMRGridPatch


class StretchedGrid(AMRGridPatch):
    def __init__(self, id, cell_widths, *, filename=None, index=None):
        self.cell_widths = [np.array(_) for _ in cell_widths]
        super().__init__(id, filename, index)

    def _check_consistency(self):
        computed_right_edge = self.LeftEdge + [
            _.sum() for _ in self.cell_widths * self.LeftEdge.uq
        ]
        assert (computed_right_edge == self.RightEdge).all()

    def _get_selector_mask(self, selector):
        if self._cache_mask and hash(selector) == self._last_selector_id:
            mask = self._last_mask
        else:
            mask = selector.fill_mask(self)
            if self._cache_mask:
                self._last_mask = mask
            self._last_selector_id = hash(selector)
            if mask is None:
                self._last_count = 0
            else:
                self._last_count = mask.sum()
        return mask

    def select_fwidth(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty((0, 3), dtype="float64")
        indices = convert_mask_to_indices(mask, self._last_count)

        coords = np.array(
            [
                self.cell_widths[0][indices[:, 0]],
                self.cell_widths[1][indices[:, 1]],
                self.cell_widths[2][indices[:, 2]],
            ]
        ).T
        return coords

    def select_fcoords(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty((0, 3), dtype="float64")
        cell_centers = [
            self.LeftEdge[i].d
            + np.add.accumulate(self.cell_widths[i])
            - 0.5 * self.cell_widths[i]
            for i in range(3)
        ]
        indices = convert_mask_to_indices(mask, self._last_count)
        coords = np.array(
            [
                cell_centers[0][indices[:, 0]],
                cell_centers[1][indices[:, 1]],
                cell_centers[2][indices[:, 2]],
            ]
        ).T
        return coords

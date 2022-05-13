import numpy as np


class FluidOperator:
    def apply(self, ds):
        for g in ds.index.grids:
            self(g)


class TopHatSphere(FluidOperator):
    def __init__(self, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields

    def __call__(self, grid, sub_select=None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax].to_ndarray() - self.center[i]) ** 2.0, r)
        np.sqrt(r, r)
        ind = r <= self.radius
        if sub_select is not None:
            ind &= sub_select
        for field, val in self.fields.items():
            grid[field][ind] = val


class CoredSphere(FluidOperator):
    def __init__(self, core_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.core_radius = core_radius

    def __call__(self, grid, sub_select=None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        r2 = self.radius**2
        cr2 = self.core_radius**2
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax] - self.center[i]) ** 2.0, r)
        np.maximum(r, cr2, r)
        ind = r <= r2
        if sub_select is not None:
            ind &= sub_select
        for field, (outer_val, inner_val) in self.fields.items():
            val = ((r[ind] - cr2) / (r2 - cr2)) ** 0.5 * (outer_val - inner_val)
            grid[field][ind] = val + inner_val


class BetaModelSphere(FluidOperator):
    def __init__(self, beta, core_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.core_radius = core_radius
        self.beta = beta

    def __call__(self, grid, sub_select=None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        r2 = self.radius**2
        cr2 = self.core_radius**2
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax].ndarray_view() - self.center[i]) ** 2.0, r)
        ind = r <= r2
        if sub_select is not None:
            ind &= sub_select
        for field, core_val in self.fields.items():
            val = core_val * (1.0 + r[ind] / cr2) ** (-1.5 * self.beta)
            grid[field][ind] = val


class NFWModelSphere(FluidOperator):
    def __init__(self, scale_radius, radius, center, fields):
        self.radius = radius
        self.center = center
        self.fields = fields
        self.scale_radius = scale_radius  # unitless

    def __call__(self, grid, sub_select=None):
        r = np.zeros(grid.ActiveDimensions, dtype="float64")
        for i, ax in enumerate("xyz"):
            np.add(r, (grid[ax].ndarray_view() - self.center[i]) ** 2.0, r)
        np.sqrt(r, r)
        ind = r <= self.radius
        r /= self.scale_radius
        if sub_select is not None:
            ind &= sub_select
        for field, scale_val in self.fields.items():
            val = scale_val / (r[ind] * (1.0 + r[ind]) ** 2)
            grid[field][ind] = val


class RandomFluctuation(FluidOperator):
    def __init__(self, fields):
        self.fields = fields

    def __call__(self, grid, sub_select=None):
        if sub_select is None:
            sub_select = Ellipsis
        for field, mag in self.fields.items():
            vals = grid[field][sub_select]
            rc = 1.0 + (np.random.random(vals.shape) - 0.5) * mag
            grid[field][sub_select] *= rc

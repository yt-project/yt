"""
Quantities that can be derived from Enzo data that may also required additional
arguments.  (Standard arguments -- such as the center of a distribution of
points -- are excluded here, and left to the EnzoDerivedFields.)



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import *

from yt.config import ytcfg
from yt.units.yt_array import YTArray, uconcatenate, array_like_field
from yt.utilities.exceptions import YTFieldNotFound
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_objects
from yt.utilities.lib.Octree import Octree
from yt.utilities.physical_constants import \
    gravitational_constant_cgs, \
    HUGE
from yt.utilities.math_utils import prec_accum

derived_quantity_registry = {}

class DerivedQuantity(ParallelAnalysisInterface):
    num_vals = -1

    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if name != "DerivedQuantity":
                derived_quantity_registry[name] = cls

    def __init__(self, data_source):
        self.data_source = data_source

    def count_values(self, *args, **kwargs):
        return

    def __call__(self, *args, **kwargs):
        """Calculate results for the derived quantity"""
        self.count_values(*args, **kwargs)
        chunks = self.data_source.chunks([], chunking_style="io")
        storage = {}
        for sto, ds in parallel_objects(chunks, -1, storage = storage):
            sto.result = self.process_chunk(ds, *args, **kwargs)
        # Now storage will have everything, and will be done via pickling, so
        # the units will be preserved.  (Credit to Nathan for this
        # idea/implementation.)
        values = [ [] for i in range(self.num_vals) ]
        for key in sorted(storage):
            for i in range(self.num_vals):
                values[i].append(storage[key][i])
        # These will be YTArrays
        values = [self.data_source.ds.arr(values[i]) for i in range(self.num_vals)]
        values = self.reduce_intermediate(values)
        return values

    def process_chunk(self, data, *args, **kwargs):
        raise NotImplementedError

    def reduce_intermediate(self, values):
        raise NotImplementedError

class DerivedQuantityCollection(object):
    def __new__(cls, data_source, *args, **kwargs):
        inst = object.__new__(cls)
        inst.data_source = data_source
        for f in inst.keys():
            setattr(inst, camelcase_to_underscore(f), inst[f])
        return inst

    def __getitem__(self, key):
        dq = derived_quantity_registry[key]
        # Instantiate here, so we can pass it the data object
        # Note that this means we instantiate every time we run help, etc
        # I have made my peace with this.
        return dq(self.data_source)

    def keys(self):
        return derived_quantity_registry.keys()

class WeightedAverageQuantity(DerivedQuantity):
    r"""
    Calculates the weight average of a field or fields.

    Where f is the field and w is the weight, the weighted average is
    Sum_i(f_i \* w_i) / Sum_i(w_i).

    Parameters
    ----------
    fields : field or list of fields
        The field or fields of which the average value is to be calculated.
    weight : field
        The weight field.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.weighted_average_quantity([("gas", "density"),
    ...                                                ("gas", "temperature")],
    ...                                               ("gas", "cell_mass"))

    """
    def count_values(self, fields, weight):
        # This is a list now
        self.num_vals = len(fields) + 1

    def __call__(self, fields, weight):
        fields = ensure_list(fields)
        rv = super(WeightedAverageQuantity, self).__call__(fields, weight)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, fields, weight):
        vals = [(data[field] * data[weight]).sum(dtype=np.float64)
                for field in fields]
        wv = data[weight].sum(dtype=np.float64)
        return vals + [wv]

    def reduce_intermediate(self, values):
        w = values.pop(-1).sum(dtype=np.float64)
        return [v.sum(dtype=np.float64)/w for v in values]

class TotalQuantity(DerivedQuantity):
    r"""
    Calculates the sum of the field or fields.

    Parameters
    ----------
    fields : field or list of fields
        The field to be summed.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.total_quantity([("gas", "cell_mass")])

    """
    def count_values(self, fields):
        # This is a list now
        self.num_vals = len(fields)

    def __call__(self, fields):
        fields = ensure_list(fields)
        rv = super(TotalQuantity, self).__call__(fields)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, fields):
        vals = [data[field].sum(dtype=np.float64)
                for field in fields]
        return vals

    def reduce_intermediate(self, values):
        return [v.sum(dtype=np.float64) for v in values]

class TotalMass(TotalQuantity):
    r"""
    Calculates the total mass in gas and particles. Returns a tuple where the
    first part is total gas mass and the second part is total particle mass.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.total_mass()

    """
    def __call__(self):
        self.data_source.ds.index
        fi = self.data_source.ds.field_info
        fields = []
        if ("gas", "cell_mass") in fi:
            fields.append(("gas", "cell_mass"))
        if ("all", "particle_mass") in fi:
            fields.append(("all", "particle_mass"))
        rv = super(TotalMass, self).__call__(fields)
        return rv

class CenterOfMass(DerivedQuantity):
    r"""
    Calculates the center of mass, using gas and/or particles.

    The center of mass is the mass-weighted mean position.

    Parameters
    ----------
    use_gas : bool
        Flag to include gas in the calculation.  Gas is ignored if not
        present.
        Default: True
    use_particles : bool
        Flag to include particles in the calculation.  Particles are ignored
        if not present.
        Default: False

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.center_of_mass()

    """
    def count_values(self, use_gas = True, use_particles = False):
        use_gas &= \
          (("gas", "cell_mass") in self.data_source.ds.field_info)
        use_particles &= \
          (("all", "particle_mass") in self.data_source.ds.field_info)
        self.num_vals = 0
        if use_gas:
            self.num_vals += 4
        if use_particles:
            self.num_vals += 4

    def process_chunk(self, data, use_gas = True, use_particles = False):
        use_gas &= \
          (("gas", "cell_mass") in self.data_source.ds.field_info)
        use_particles &= \
          (("all", "particle_mass") in self.data_source.ds.field_info)
        vals = []
        if use_gas:
            vals += [(data[ax] * data["gas", "cell_mass"]).sum(dtype=np.float64)
                     for ax in 'xyz']
            vals.append(data["gas", "cell_mass"].sum(dtype=np.float64))
        if use_particles:
            vals += [(data["all", "particle_position_%s" % ax] *
                      data["all", "particle_mass"]).sum(dtype=np.float64)
                     for ax in 'xyz']
            vals.append(data["all", "particle_mass"].sum(dtype=np.float64))
        return vals

    def reduce_intermediate(self, values):
        if len(values) not in (4, 8):
            raise RuntimeError
        x = values.pop(0).sum(dtype=np.float64)
        y = values.pop(0).sum(dtype=np.float64)
        z = values.pop(0).sum(dtype=np.float64)
        w = values.pop(0).sum(dtype=np.float64)
        if len(values) > 0:
            # Note that this could be shorter if we pre-initialized our x,y,z,w
            # values as YTQuantity objects.
            x += values.pop(0).sum(dtype=np.float64)
            y += values.pop(0).sum(dtype=np.float64)
            z += values.pop(0).sum(dtype=np.float64)
            w += values.pop(0).sum(dtype=np.float64)
        return self.data_source.ds.arr([v/w for v in [x, y, z]])

class BulkVelocity(DerivedQuantity):
    r"""
    Calculates the bulk velocity, using gas and/or particles.

    The bulk velocity is the mass-weighted mean velocity.

    Parameters
    ----------
    use_gas : bool
        Flag to include gas in the calculation.  Gas is ignored if not
        present.
        Default: True
    use_particles : bool
        Flag to include particles in the calculation.  Particles are ignored
        if not present.
        Default: True

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.bulk_velocity()

    """
    def count_values(self, use_gas = True, use_particles = False):
        # This is a list now
        self.num_vals = 0
        if use_gas:
            self.num_vals += 4
        if use_particles:
            self.num_vals += 4

    def process_chunk(self, data, use_gas = True, use_particles = False):
        vals = []
        if use_gas:
            vals += [(data["gas", "velocity_%s" % ax] *
                      data["gas", "cell_mass"]).sum(dtype=np.float64)
                     for ax in 'xyz']
            vals.append(data["gas", "cell_mass"].sum(dtype=np.float64))
        if use_particles:
            vals += [(data["all", "particle_velocity_%s" % ax] *
                      data["all", "particle_mass"]).sum(dtype=np.float64)
                     for ax in 'xyz']
            vals.append(data["all", "particle_mass"].sum(dtype=np.float64))
        return vals

    def reduce_intermediate(self, values):
        if len(values) not in (4, 8):
            raise RuntimeError
        x = values.pop(0).sum(dtype=np.float64)
        y = values.pop(0).sum(dtype=np.float64)
        z = values.pop(0).sum(dtype=np.float64)
        w = values.pop(0).sum(dtype=np.float64)
        if len(values) > 0:
            # Note that this could be shorter if we pre-initialized our x,y,z,w
            # values as YTQuantity objects.
            x += values.pop(0).sum(dtype=np.float64)
            y += values.pop(0).sum(dtype=np.float64)
            z += values.pop(0).sum(dtype=np.float64)
            w += values.pop(0).sum(dtype=np.float64)
        return self.data_source.ds.arr([v/w for v in [x, y, z]])

class WeightedVariance(DerivedQuantity):
    r"""
    Calculates the weighted variance and weighted mean for a field
    or list of fields.

    Where f is the field, w is the weight, and <f_w> is the weighted mean,
    the weighted variance is
    Sum_i( (f_i - <f_w>)^2 \* w_i ) / Sum_i(w_i).

    Parameters
    ----------
    fields : field or list of fields
        The field or fields of which the variance and mean values are
        to be calculated.
    weight : field
        The weight field.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.weighted_variance([("gas", "density"),
    ...                                        ("gas", "temperature")],
    ...                                       ("gas", "cell_mass"))

    """
    def count_values(self, fields, weight):
        # This is a list now
        self.num_vals = 2 * len(fields) + 1

    def __call__(self, fields, weight):
        fields = ensure_list(fields)
        rv = super(WeightedVariance, self).__call__(fields, weight)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, fields, weight):
        my_weight = data[weight].sum(dtype=np.float64)
        if my_weight == 0:
            return [0.0 for field in fields] + \
              [0.0 for field in fields] + [0.0]
        my_means = [(data[field] *  data[weight]).sum(dtype=np.float64) / my_weight
                    for field in fields]
        my_var2s = [(data[weight] * (data[field] -
                                     my_mean)**2).sum(dtype=np.float64) / my_weight
                   for field, my_mean in zip(fields, my_means)]
        return my_means + my_var2s + [my_weight]

    def reduce_intermediate(self, values):
        my_weight = values.pop(-1)
        all_weight = my_weight.sum(dtype=np.float64)
        rvals = []
        for i in range(len(values) / 2):
            my_mean = values[i]
            my_var2 = values[i + len(values) / 2]
            all_mean = (my_weight * my_mean).sum(dtype=np.float64) / all_weight
            rvals.append(np.sqrt((my_weight * (my_var2 +
                                               (my_mean - all_mean)**2)).sum(dtype=np.float64) /
                                               all_weight))
            rvals.append(all_mean)
        return rvals

class AngularMomentumVector(DerivedQuantity):
    r"""
    Calculates the angular momentum vector, using gas and/or particles.

    The angular momentum vector is the mass-weighted mean specific angular momentum.

    Parameters
    ----------
    use_gas : bool
        Flag to include gas in the calculation.  Gas is ignored if not
        present.
        Default: True
    use_particles : bool
        Flag to include particles in the calculation.  Particles are ignored
        if not present.
        Default: True

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.angular_momentum_vector()

    """
    def count_values(self, use_gas=True, use_particles=True):
        use_gas &= \
          (("gas", "cell_mass") in self.data_source.ds.field_info)
        use_particles &= \
          (("all", "particle_mass") in self.data_source.ds.field_info)
        num_vals = 0
        if use_gas: num_vals += 4
        if use_particles: num_vals += 4
        self.num_vals = num_vals

    def process_chunk(self, data, use_gas=True, use_particles=True):
        use_gas &= \
          (("gas", "cell_mass") in self.data_source.ds.field_info)
        use_particles &= \
          (("all", "particle_mass") in self.data_source.ds.field_info)
        rvals = []
        if use_gas:
            rvals.extend([(data["gas", "specific_angular_momentum_%s" % axis] *
                           data["gas", "cell_mass"]).sum(dtype=np.float64) \
                          for axis in "xyz"])
            rvals.append(data["gas", "cell_mass"].sum(dtype=np.float64))
        if use_particles:
            rvals.extend([(data["all", "particle_specific_angular_momentum_%s" % axis] *
                           data["all", "particle_mass"]).sum(dtype=np.float64) \
                          for axis in "xyz"])
            rvals.append(data["all", "particle_mass"].sum(dtype=np.float64))
        return rvals

    def reduce_intermediate(self, values):
        jx = values.pop(0).sum(dtype=np.float64)
        jy = values.pop(0).sum(dtype=np.float64)
        jz = values.pop(0).sum(dtype=np.float64)
        m  = values.pop(0).sum(dtype=np.float64)
        if values:
            jx += values.pop(0).sum(dtype=np.float64)
            jy += values.pop(0).sum(dtype=np.float64)
            jz += values.pop(0).sum(dtype=np.float64)
            m  += values.pop(0).sum(dtype=np.float64)
        return (jx / m, jy / m, jz / m)

class Extrema(DerivedQuantity):
    r"""
    Calculates the min and max value of a field or list of fields.

    Parameters
    ----------
    fields : field or list of fields
        The field over which the extrema are to be calculated.
    non_zero : bool
        If True, only positive values are considered in the calculation.
        Default: False

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.extrema([("gas", "density"),
    ...                              ("gas", "temperature")])

    """
    def count_values(self, fields, non_zero):
        self.num_vals = len(fields) * 2

    def __call__(self, fields, non_zero = False):
        fields = ensure_list(fields)
        rv = super(Extrema, self).__call__(fields, non_zero)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, fields, non_zero):
        vals = []
        for field in fields:
            field = data._determine_fields(field)[0]
            fd = data[field]
            if non_zero: fd = fd[fd > 0.0]
            if fd.size > 0:
                vals += [fd.min(), fd.max()]
            else:
                vals += [array_like_field(data, HUGE, field),
                         array_like_field(data, -HUGE, field)]
        return vals

    def reduce_intermediate(self, values):
        # The values get turned into arrays here.
        return [(mis.min(), mas.max() )
                for mis, mas in zip(values[::2], values[1::2])]

class MaxLocation(DerivedQuantity):
    r"""
    Calculates the maximum value plus the index, x, y, and z position
    of the maximum.

    Parameters
    ----------
    field : field
        The field over which the extrema are to be calculated.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.max_location(("gas", "density"))

    """
    def count_values(self, *args, **kwargs):
        self.num_vals = 5

    def __call__(self, field):
        rv = super(MaxLocation, self).__call__(field)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, field):
        field = data._determine_fields(field)[0]
        ma = array_like_field(data, -HUGE, field)
        mx = array_like_field(data, -1, "x")
        my = array_like_field(data, -1, "y")
        mz = array_like_field(data, -1, "z")
        maxi = -1
        if data[field].size > 0:
            maxi = np.argmax(data[field])
            ma = data[field][maxi]
            mx, my, mz = [data[ax][maxi] for ax in 'xyz']
        return (ma, maxi, mx, my, mz)

    def reduce_intermediate(self, values):
        i = np.argmax(values[0]) # ma is values[0]
        return [val[i] for val in values]

class MinLocation(DerivedQuantity):
    r"""
    Calculates the minimum value plus the index, x, y, and z position
    of the minimum.

    Parameters
    ----------
    field : field
        The field over which the extrema are to be calculated.

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.min_location(("gas", "density"))

    """
    def count_values(self, *args, **kwargs):
        self.num_vals = 5

    def __call__(self, field):
        rv = super(MinLocation, self).__call__(field)
        if len(rv) == 1: rv = rv[0]
        return rv

    def process_chunk(self, data, field):
        field = data._determine_fields(field)[0]
        ma = array_like_field(data, HUGE, field)
        mx = array_like_field(data, -1, "x")
        my = array_like_field(data, -1, "y")
        mz = array_like_field(data, -1, "z")
        mini = -1
        if data[field].size > 0:
            mini = np.argmin(data[field])
            ma = data[field][mini]
            mx, my, mz = [data[ax][mini] for ax in 'xyz']
        return (ma, mini, mx, my, mz)

    def reduce_intermediate(self, values):
        i = np.argmin(values[0]) # ma is values[0]
        return [val[i] for val in values]

class SpinParameter(DerivedQuantity):
    r"""
    Calculates the dimensionless spin parameter.

    Given by Equation 3 of Peebles (1971, A&A, 11, 377), the spin parameter
    is defined as

    lambda = (L \* |E|^(1/2)) / (G \* M^5/2),

    where L is the total angular momentum, E is the total energy (kinetic and
    potential), G is the gravitational constant, and M is the total mass.

    Parameters
    ----------
    use_gas : bool
        Flag to include gas in the calculation.  Gas is ignored if not
        present.
        Default: True
    use_particles : bool
        Flag to include particles in the calculation.  Particles are ignored
        if not present.
        Default: True

    Examples
    --------

    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> print ad.quantities.center_of_mass()

    """
    def count_values(self, **kwargs):
        self.num_vals = 3

    def process_chunk(self, data, use_gas=True, use_particles=True):
        use_gas &= \
          (("gas", "cell_mass") in self.data_source.ds.field_info)
        use_particles &= \
          (("all", "particle_mass") in self.data_source.ds.field_info)
        e = data.ds.quan(0., "erg")
        j = data.ds.quan(0., "g*cm**2/s")
        m = data.ds.quan(0., "g")
        if use_gas:
            e += (data["gas", "kinetic_energy"] *
                  data["index", "cell_volume"]).sum(dtype=np.float64)
            j += data["gas", "angular_momentum_magnitude"].sum(dtype=np.float64)
            m += data["gas", "cell_mass"].sum(dtype=np.float64)
        if use_particles:
            e += (data["all", "particle_velocity_magnitude"]**2 *
                  data["all", "particle_mass"]).sum(dtype=np.float64)
            j += data["all", "particle_angular_momentum_magnitude"].sum(dtype=np.float64)
            m += data["all", "particle_mass"].sum(dtype=np.float64)
        return (e, j, m)

    def reduce_intermediate(self, values):
        e = values.pop(0).sum(dtype=np.float64)
        j = values.pop(0).sum(dtype=np.float64)
        m = values.pop(0).sum(dtype=np.float64)
        return j * np.sqrt(np.abs(e)) / m**2.5 / gravitational_constant_cgs

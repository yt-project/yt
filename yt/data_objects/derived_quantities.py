"""
Quantities that can be derived from Enzo data that may also required additional
arguments.  (Standard arguments -- such as the center of a distribution of
points -- are excluded here, and left to the EnzoDerivedFields.)

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as na

from yt.funcs import *

from yt.config import ytcfg
from yt.data_objects.field_info_container import \
    FieldDetector
from yt.utilities.data_point_utilities import FindBindingEnergy
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.lib import Octree

__CUDA_BLOCK_SIZE = 256

quantity_info = {}

class GridChildMaskWrapper:
    def __init__(self, grid, data_source):
        self.grid = grid
        self.data_source = data_source
        # We have a local cache so that *within* a call to the DerivedQuantity
        # function, we only read each field once.  Otherwise, when preloading
        # the field will be popped and removed and lost if the underlying data
        # source's _get_data_from_grid method is wrapped by restore_state.
        # This is common.  So, if data[something] is accessed multiple times by
        # a single quantity, the second time will re-read the data the slow
        # way.
        self.local_cache = {}
    def __getattr__(self, attr):
        return getattr(self.grid, attr)
    def __getitem__(self, item):
        if item not in self.local_cache:
            data = self.data_source._get_data_from_grid(self.grid, item)
            self.local_cache[item] = data
        return self.local_cache[item]

class DerivedQuantity(ParallelAnalysisInterface):
    def __init__(self, collection, name, function,
                 combine_function, units = "",
                 n_ret = 0, force_unlazy=False):
        # We wrap the function with our object
        ParallelAnalysisInterface.__init__(self)
        self.__doc__ = function.__doc__
        self.__name__ = name
        self.collection = collection
        self._data_source = collection.data_source
        self.func = function
        self.c_func = combine_function
        self.n_ret = n_ret
        self.force_unlazy = force_unlazy

    def __call__(self, *args, **kwargs):
        lazy_reader = kwargs.pop('lazy_reader', True)
        preload = kwargs.pop('preload', ytcfg.getboolean("yt","__parallel"))
        if preload:
            if not lazy_reader: mylog.debug("Turning on lazy_reader because of preload")
            lazy_reader = True
            e = FieldDetector(flat = True)
            e.NumberOfParticles = 1
            self.func(e, *args, **kwargs)
            mylog.debug("Preloading %s", e.requested)
            self.comm.preload([g for g in self._get_grid_objs()], e.requested,
                          self._data_source.pf.h.io)
        if lazy_reader and not self.force_unlazy:
            return self._call_func_lazy(args, kwargs)
        else:
            return self._call_func_unlazy(args, kwargs)

    def _call_func_lazy(self, args, kwargs):
        self.retvals = [ [] for i in range(self.n_ret)]
        for gi,g in enumerate(self._get_grids()):
            rv = self.func(GridChildMaskWrapper(g, self._data_source), *args, **kwargs)
            if not iterable(rv): rv = (rv,)
            for i in range(self.n_ret): self.retvals[i].append(rv[i])
            g.clear_data()
        self.retvals = [na.array(self.retvals[i]) for i in range(self.n_ret)]
        return self.c_func(self._data_source, *self.retvals)

    def _finalize_parallel(self):
        # Note that we do some fancy footwork here.
        # _par_combine_object and its affiliated alltoall function
        # assume that the *long* axis is the last one.  However,
        # our long axis is the first one!
        rv = []
        for my_list in self.retvals:
            data = na.array(my_list).transpose()
            rv.append(self.comm.par_combine_object(data,
                        datatype="array", op="cat").transpose())
        self.retvals = rv
        
    def _call_func_unlazy(self, args, kwargs):
        retval = self.func(self._data_source, *args, **kwargs)
        return self.c_func(self._data_source, *retval)

def add_quantity(name, **kwargs):
    if 'function' not in kwargs or 'combine_function' not in kwargs:
        mylog.error("Not adding field %s because both function and combine_function must be provided" % name)
        return
    f = kwargs.pop('function')
    c = kwargs.pop('combine_function')
    quantity_info[name] = (name, f, c, kwargs)

class DerivedQuantityCollection(object):
    functions = quantity_info
    def __init__(self, data_source):
        self.data_source = data_source

    def __getitem__(self, key):
        if key not in self.functions:
            raise KeyError(key)
        args = self.functions[key][:3]
        kwargs = self.functions[key][3]
        # Instantiate here, so we can pass it the data object
        # Note that this means we instantiate every time we run help, etc
        # I have made my peace with this.
        return DerivedQuantity(self, *args, **kwargs)

    def keys(self):
        return self.functions.keys()

def _TotalMass(data):
    """
    This function takes no arguments and returns the sum of cell masses and
    particle masses in the object.
    """
    baryon_mass = data["CellMassMsun"].sum()
    particle_mass = data["ParticleMassMsun"].sum()
    return [baryon_mass + particle_mass]
def _combTotalMass(data, total_mass):
    return total_mass.sum()
add_quantity("TotalMass", function=_TotalMass,
             combine_function=_combTotalMass, n_ret=1)

def _CenterOfMass(data, use_cells=True, use_particles=False):
    """
    This function returns the location of the center
    of mass. By default, it computes of the *non-particle* data in the object. 

    Parameters
    ----------

    use_cells : bool
        If True, will include the cell mass (default: True)
    use_particles : bool
        if True, will include the particles in the object (default: False)
    """
    x = y = z = den = 0
    if use_cells: 
        x += (data["x"] * data["CellMassMsun"]).sum()
        y += (data["y"] * data["CellMassMsun"]).sum()
        z += (data["z"] * data["CellMassMsun"]).sum()
        den += data["CellMassMsun"].sum()
    if use_particles:
        x += (data["particle_position_x"] * data["ParticleMassMsun"]).sum()
        y += (data["particle_position_y"] * data["ParticleMassMsun"]).sum()
        z += (data["particle_position_z"] * data["ParticleMassMsun"]).sum()
        den += data["ParticleMassMsun"].sum()

    return x,y,z, den
def _combCenterOfMass(data, x,y,z, den):
    return na.array([x.sum(), y.sum(), z.sum()])/den.sum()
add_quantity("CenterOfMass", function=_CenterOfMass,
             combine_function=_combCenterOfMass, n_ret = 4)

def _WeightedAverageQuantity(data, field, weight):
    """
    This function returns an averaged quantity.

    :param field: The field to average
    :param weight: The field to weight by
    """
    num = (data[field] * data[weight]).sum()
    den = data[weight].sum()
    return num, den
def _combWeightedAverageQuantity(data, field, weight):
    return field.sum()/weight.sum()
add_quantity("WeightedAverageQuantity", function=_WeightedAverageQuantity,
             combine_function=_combWeightedAverageQuantity, n_ret = 2)

def _BulkVelocity(data):
    """
    This function returns the mass-weighted average velocity in the object.
    """
    xv = (data["x-velocity"] * data["CellMassMsun"]).sum()
    yv = (data["y-velocity"] * data["CellMassMsun"]).sum()
    zv = (data["z-velocity"] * data["CellMassMsun"]).sum()
    w = data["CellMassMsun"].sum()
    return xv, yv, zv, w
def _combBulkVelocity(data, xv, yv, zv, w):
    w = w.sum()
    xv = xv.sum()/w
    yv = yv.sum()/w
    zv = zv.sum()/w
    return na.array([xv, yv, zv])
add_quantity("BulkVelocity", function=_BulkVelocity,
             combine_function=_combBulkVelocity, n_ret=4)

def _AngularMomentumVector(data):
    """
    This function returns the mass-weighted average angular momentum vector.
    """
    amx = data["SpecificAngularMomentumX"]*data["CellMassMsun"]
    amy = data["SpecificAngularMomentumY"]*data["CellMassMsun"]
    amz = data["SpecificAngularMomentumZ"]*data["CellMassMsun"]
    j_mag = [amx.sum(), amy.sum(), amz.sum()]
    return [j_mag]

def _StarAngularMomentumVector(data):
    """
    This function returns the mass-weighted average angular momentum vector 
    for stars.
    """
    is_star = data["creation_time"] > 0
    star_mass = data["ParticleMassMsun"][is_star]
    sLx = data["ParticleSpecificAngularMomentumX"][is_star]
    sLy = data["ParticleSpecificAngularMomentumY"][is_star]
    sLz = data["ParticleSpecificAngularMomentumZ"][is_star]
    amx = sLx * star_mass
    amy = sLy * star_mass
    amz = sLz * star_mass
    j_mag = [amx.sum(), amy.sum(), amz.sum()]
    return [j_mag]

def _combAngularMomentumVector(data, j_mag):
    if len(j_mag.shape) < 2: j_mag = na.expand_dims(j_mag, 0)
    L_vec = j_mag.sum(axis=0)
    L_vec_norm = L_vec / na.sqrt((L_vec**2.0).sum())
    return L_vec_norm
add_quantity("AngularMomentumVector", function=_AngularMomentumVector,
             combine_function=_combAngularMomentumVector, n_ret=1)

add_quantity("StarAngularMomentumVector", function=_StarAngularMomentumVector,
             combine_function=_combAngularMomentumVector, n_ret=1)

def _BaryonSpinParameter(data):
    """
    This function returns the spin parameter for the baryons, but it uses
    the particles in calculating enclosed mass.
    """
    m_enc = data["CellMassMsun"].sum() + data["ParticleMassMsun"].sum()
    amx = data["SpecificAngularMomentumX"]*data["CellMassMsun"]
    amy = data["SpecificAngularMomentumY"]*data["CellMassMsun"]
    amz = data["SpecificAngularMomentumZ"]*data["CellMassMsun"]
    j_mag = na.array([amx.sum(), amy.sum(), amz.sum()])
    e_term_pre = na.sum(data["CellMassMsun"]*data["VelocityMagnitude"]**2.0)
    weight=data["CellMassMsun"].sum()
    return j_mag, m_enc, e_term_pre, weight
def _combBaryonSpinParameter(data, j_mag, m_enc, e_term_pre, weight):
    # Because it's a vector field, we have to ensure we have enough dimensions
    if len(j_mag.shape) < 2: j_mag = na.expand_dims(j_mag, 0)
    W = weight.sum()
    M = m_enc.sum()
    J = na.sqrt(((j_mag.sum(axis=0))**2.0).sum())/W
    E = na.sqrt(e_term_pre.sum()/W)
    G = 6.67e-8 # cm^3 g^-1 s^-2
    spin = J * E / (M*1.989e33*G)
    return spin
add_quantity("BaryonSpinParameter", function=_BaryonSpinParameter,
             combine_function=_combBaryonSpinParameter, n_ret=4)

def _ParticleSpinParameter(data):
    """
    This function returns the spin parameter for the baryons, but it uses
    the particles in calculating enclosed mass.
    """
    m_enc = data["CellMassMsun"].sum() + data["ParticleMassMsun"].sum()
    amx = data["ParticleSpecificAngularMomentumX"]*data["ParticleMassMsun"]
    if amx.size == 0: return (na.zeros((3,), dtype='float64'), m_enc, 0, 0)
    amy = data["ParticleSpecificAngularMomentumY"]*data["ParticleMassMsun"]
    amz = data["ParticleSpecificAngularMomentumZ"]*data["ParticleMassMsun"]
    j_mag = na.array([amx.sum(), amy.sum(), amz.sum()])
    e_term_pre = na.sum(data["ParticleMassMsun"]
                       *data["ParticleVelocityMagnitude"]**2.0)
    weight=data["ParticleMassMsun"].sum()
    return j_mag, m_enc, e_term_pre, weight
add_quantity("ParticleSpinParameter", function=_ParticleSpinParameter,
             combine_function=_combBaryonSpinParameter, n_ret=4)
    
def _IsBound(data, truncate = True, include_thermal_energy = False,
    treecode = True, opening_angle = 1.0, periodic_test = False, include_particles = True):
    r"""
    This returns whether or not the object is gravitationally bound. If this
    returns a value greater than one, it is bound, and otherwise not.
    
    Parameters
    ----------
    truncate : Bool
        Should the calculation stop once the ratio of
        gravitational:kinetic is 1.0?
    include_thermal_energy : Bool
        Should we add the energy from ThermalEnergy
        on to the kinetic energy to calculate 
        binding energy?
    treecode : Bool
        Whether or not to use the treecode.
    opening_angle : Float 
        The maximal angle a remote node may subtend in order
        for the treecode method of mass conglomeration may be
        used to calculate the potential between masses.
    periodic_test : Bool 
        Used for testing the periodic adjustment machinery
        of this derived quantity.
    include_particles : Bool
	Should we add the mass contribution of particles
	to calculate binding energy?

    Examples
    --------
    >>> sp.quantities["IsBound"](truncate=False,
    ... include_thermal_energy=True, treecode=False, opening_angle=2.0)
    0.32493
    """
    # Kinetic energy
    bv_x,bv_y,bv_z = data.quantities["BulkVelocity"]()
    # One-cell objects are NOT BOUND.
    if data["CellMass"].size == 1: return [0.0]

    kinetic = 0.5 * (data["CellMass"] * 
                     ((data["x-velocity"] - bv_x)**2 + 
                      (data["y-velocity"] - bv_y)**2 +
                      (data["z-velocity"] - bv_z)**2)).sum()

    if (include_particles):
	mass_to_use = data["TotalMass"]
        kinetic += 0.5 * (data["Dark_Matter_Mass"] *
                          ((data["cic_particle_velocity_x"] - bv_x)**2 +
                           (data["cic_particle_velocity_y"] - bv_y)**2 +
                           (data["cic_particle_velocity_z"] - bv_z)**2)).sum()
    else:
	mass_to_use = data["CellMass"]
    # Add thermal energy to kinetic energy
    if (include_thermal_energy):
        thermal = (data["ThermalEnergy"] * mass_to_use).sum()
        kinetic += thermal
    if periodic_test:
        kinetic = na.ones_like(kinetic)
    # Gravitational potential energy
    # We only divide once here because we have velocity in cgs, but radius is
    # in code.
    G = 6.67e-8 / data.convert("cm") # cm^3 g^-1 s^-2
    # Check for periodicity of the clump.
    two_root = 2. / na.array(data.pf.domain_dimensions)
    domain_period = data.pf.domain_right_edge - data.pf.domain_left_edge
    periodic = na.array([0., 0., 0.])
    for i,dim in enumerate(["x", "y", "z"]):
        sorted = data[dim][data[dim].argsort()]
        # If two adjacent values are different by (more than) two root grid
        # cells, I think it's reasonable to assume that the clump wraps around.
        diff = sorted[1:] - sorted[0:-1]
        if (diff >= two_root[i]).any():
            mylog.info("Adjusting clump for periodic boundaries in dim %s" % dim)
            # We will record the distance of the larger of the two values that
            # define the gap from the right boundary, which we'll use for the
            # periodic adjustment later.
            sel = (diff >= two_root[i])
            index = na.min(na.nonzero(sel))
            # The last addition term below ensures that the data makes a full
            # wrap-around.
            periodic[i] = data.pf.domain_right_edge[i] - sorted[index + 1] + \
                two_root[i] / 2.
    # This dict won't make a copy of the data, but it will make a copy to 
    # change if needed in the periodic section immediately below.
    local_data = {}
    for label in ["x", "y", "z"]: # Separating CellMass from the for loop
        local_data[label] = data[label]
    local_data["CellMass"] = mass_to_use # Adding CellMass separately
					 # NOTE: if include_particles = True, local_data["CellMass"]
					 #       is not the same as data["CellMass"]!!!
    if periodic.any():
        # Adjust local_data to re-center the clump to remove the periodicity
        # by the gap calculated above.
        for i,dim in enumerate(["x", "y", "z"]):
            if not periodic[i]: continue
            local_data[dim] = data[dim].copy()
            local_data[dim] += periodic[i]
            local_data[dim] %= domain_period[i]
    if periodic_test:
        local_data["CellMass"] = na.ones_like(local_data["CellMass"])
    import time
    t1 = time.time()
    if treecode:
        # Calculate the binding energy using the treecode method.
        # Faster but less accurate.
        # The octree doesn't like uneven root grids, so we will make it cubical.
        root_dx = 1./na.array(data.pf.domain_dimensions).astype('float64')
        left = min([na.amin(local_data['x']), na.amin(local_data['y']),
            na.amin(local_data['z'])])
        right = max([na.amax(local_data['x']), na.amax(local_data['y']),
            na.amax(local_data['z'])])
        cover_min = na.array([left, left, left])
        cover_max = na.array([right, right, right])
        # Fix the coverage to match to root grid cell left 
        # edges for making indexes.
        cover_min = cover_min - cover_min % root_dx
        cover_max = cover_max - cover_max % root_dx
        cover_imin = (cover_min * na.array(data.pf.domain_dimensions)).astype('int64')
        cover_imax = (cover_max * na.array(data.pf.domain_dimensions) + 1).astype('int64')
        cover_ActiveDimensions = cover_imax - cover_imin
        # Create the octree with these dimensions.
        # One value (mass) with incremental=True.
        octree = Octree(cover_ActiveDimensions, 1, True)
        #print 'here', cover_ActiveDimensions
        # Now discover what levels this data comes from, not assuming
        # symmetry.
        dxes = na.unique(data['dx']) # unique returns a sorted array,
        dyes = na.unique(data['dy']) # so these will all have the same
        dzes = na.unique(data['dz']) # order.
        # We only need one dim to figure out levels, we'll use x.
        dx = 1./data.pf.domain_dimensions[0]
        levels = (na.log(dx / dxes) / na.log(data.pf.refine_by)).astype('int')
        lsort = levels.argsort()
        levels = levels[lsort]
        dxes = dxes[lsort]
        dyes = dyes[lsort]
        dzes = dzes[lsort]
        # Now we add actual data to the octree.
        for L, dx, dy, dz in zip(levels, dxes, dyes, dzes):
            mylog.info("Adding data to Octree for level %d" % L)
            sel = (data["dx"] == dx)
            thisx = (local_data["x"][sel] / dx).astype('int64') - cover_imin[0] * 2**L
            thisy = (local_data["y"][sel] / dy).astype('int64') - cover_imin[1] * 2**L
            thisz = (local_data["z"][sel] / dz).astype('int64') - cover_imin[2] * 2**L
	    vals = na.array([local_data["CellMass"][sel]], order='F')
            octree.add_array_to_tree(L, thisx, thisy, thisz, vals,
               na.ones_like(thisx).astype('float64'), treecode = 1)
        # Now we calculate the binding energy using a treecode.
        octree.finalize(treecode = 1)
        mylog.info("Using a treecode to find gravitational energy for %d cells." % local_data['x'].size)
        pot = G*octree.find_binding_energy(truncate, kinetic/G, root_dx,
            opening_angle)
        #octree.print_all_nodes()
    else:
        try:
            pot = G*_cudaIsBound(local_data, truncate, kinetic/G)
        except (ImportError, AssertionError):
            pot = G*FindBindingEnergy(local_data["CellMass"],
                                local_data['x'],local_data['y'],local_data['z'],
                                truncate, kinetic/G)
    mylog.info("Boundedness check took %0.3e seconds", time.time()-t1)
    del local_data
    return [(pot / kinetic)]
def _combIsBound(data, bound):
    return bound
add_quantity("IsBound",function=_IsBound,combine_function=_combIsBound,n_ret=1,
             force_unlazy=True)

def _cudaIsBound(data, truncate, ratio):
    bsize = __CUDA_BLOCK_SIZE
    import pycuda.driver as cuda
    import pycuda.autoinit
    import pycuda.gpuarray as gpuarray
    my_stream = cuda.Stream()
    cuda.init()
    assert cuda.Device.count() >= 1

    mass_scale_factor = 1.0/(data['CellMass'].max())
    m = (data['CellMass'] * mass_scale_factor).astype('float32')
    assert(m.size > bsize)

    gsize=int(na.ceil(float(m.size)/bsize))
    assert(gsize > 16)

    # Now the tedious process of rescaling our values...
    length_scale_factor = data['dx'].max()/data['dx'].min()
    x = ((data['x'] - data['x'].min()) * length_scale_factor).astype('float32')
    y = ((data['y'] - data['y'].min()) * length_scale_factor).astype('float32')
    z = ((data['z'] - data['z'].min()) * length_scale_factor).astype('float32')
    p = na.zeros(z.shape, dtype='float32')
    
    x_gpu = cuda.mem_alloc(x.size * x.dtype.itemsize)
    y_gpu = cuda.mem_alloc(y.size * y.dtype.itemsize)
    z_gpu = cuda.mem_alloc(z.size * z.dtype.itemsize)
    m_gpu = cuda.mem_alloc(m.size * m.dtype.itemsize)
    p_gpu = cuda.mem_alloc(p.size * p.dtype.itemsize)
    for ag, a in [(x_gpu, x), (y_gpu, y), (z_gpu, z), (m_gpu, m), (p_gpu, p)]:
        cuda.memcpy_htod(ag, a)
    source = """

      extern __shared__ float array[];

      __global__ void isbound(float *x, float *y, float *z, float *m,
                              float *p, int *nelem)
      {

        /* My index in the array */
        int idx1 = blockIdx.x * blockDim.x + threadIdx.x;
        /* Note we are setting a start index */
        int idx2 = blockIdx.y * blockDim.x;
        int offset = threadIdx.x;

        /* Here we're just setting up convenience pointers to our
           shared array */

        float* x_data1 = (float*) array;
        float* y_data1 = (float*) &x_data1[blockDim.x];
        float* z_data1 = (float*) &y_data1[blockDim.x];
        float* m_data1 = (float*) &z_data1[blockDim.x];

        float* x_data2 = (float*) &m_data1[blockDim.x];
        float* y_data2 = (float*) &x_data2[blockDim.x];
        float* z_data2 = (float*) &y_data2[blockDim.x];
        float* m_data2 = (float*) &z_data2[blockDim.x];

        x_data1[offset] = x[idx1];
        y_data1[offset] = y[idx1];
        z_data1[offset] = z[idx1];
        m_data1[offset] = m[idx1];

        x_data2[offset] = x[idx2 + offset];
        y_data2[offset] = y[idx2 + offset];
        z_data2[offset] = z[idx2 + offset];
        m_data2[offset] = m[idx2 + offset];

        __syncthreads();

        float tx, ty, tz;

        float my_p = 0.0;

        if(idx1 < %(p)s) {
            for (int i = 0; i < blockDim.x; i++){
                if(i + idx2 < idx1 + 1) continue;
                tx = (x_data1[offset]-x_data2[i]);
                ty = (y_data1[offset]-y_data2[i]);
                tz = (z_data1[offset]-z_data2[i]);
                my_p += m_data1[offset]*m_data2[i] /
                    sqrt(tx*tx+ty*ty+tz*tz);
            }
        }
        p[idx1] += my_p;
        __syncthreads();
      }
    """
    mod = cuda.SourceModule(source % dict(p=m.size))
    func = mod.get_function('isbound')
    mylog.info("Running CUDA functions.  May take a while.  (%0.5e, %s)",
               x.size, gsize)
    import pycuda.tools as ct
    t1 = time.time()
    ret = func(x_gpu, y_gpu, z_gpu, m_gpu, p_gpu,
               shared=8*bsize*m.dtype.itemsize,
         block=(bsize,1,1), grid=(gsize, gsize), time_kernel=True)
    cuda.memcpy_dtoh(p, p_gpu)
    p1 = p.sum()
    if na.any(na.isnan(p)): raise ValueError
    return p1 * (length_scale_factor / (mass_scale_factor**2.0))

def _Extrema(data, fields, non_zero = False, filter=None):
    """
    This function returns the extrema of a set of fields
    
    :param fields: A field name, or a list of field names
    :param filter: a string to be evaled to serve as a data filter.
    """
    # There is a heck of a lot of logic in this.  I really wish it were more
    # elegant.
    fields = ensure_list(fields)
    if filter is not None: this_filter = eval(filter)
    mins, maxs = [], []
    for field in fields:
        if data[field].size < 1:
            mins.append(1e90)
            maxs.append(-1e90)
            continue
        if filter is None:
            if non_zero:
                nz_filter = data[field]>0.0
                if not nz_filter.any():
                    mins.append(1e90)
                    maxs.append(-1e90)
                    continue
            else:
                nz_filter = None
            mins.append(data[field][nz_filter].min())
            maxs.append(data[field][nz_filter].max())
        else:
            if this_filter.any():
                if non_zero:
                    nz_filter = ((this_filter) &
                                 (data[field][this_filter] > 0.0))
                else: nz_filter = this_filter
                mins.append(data[field][nz_filter].min())
                maxs.append(data[field][nz_filter].max())
            else:
                mins.append(1e90)
                maxs.append(-1e90)
    return len(fields), mins, maxs
def _combExtrema(data, n_fields, mins, maxs):
    mins, maxs = na.atleast_2d(mins, maxs)
    n_fields = mins.shape[1]
    return [(na.min(mins[:,i]), na.max(maxs[:,i])) for i in range(n_fields)]
add_quantity("Extrema", function=_Extrema, combine_function=_combExtrema,
             n_ret=3)

def _Action(data, action, combine_action, filter=None):
    """
    This function evals the string given by the action arg and uses 
    the function thrown with the combine_action to combine the values.  
    A filter can be thrown to be evaled to short-circuit the calculation 
    if some criterion is not met.
    :param action: a string containing the desired action to be evaled.
    :param combine_action: the function used to combine the answers when done lazily.
    :param filter: a string to be evaled to serve as a data filter.
    """
    if filter is not None:
        if not eval(filter).any(): return 0, False, combine_action
    value = eval(action)
    return value, True, combine_action
def _combAction(data, value, valid, combine_action):
    return combine_action[0](value[valid])
add_quantity("Action", function=_Action, combine_function=_combAction, n_ret=3)

def _MaxLocation(data, field):
    """
    This function returns the location of the maximum of a set
    of fields.
    """
    ma, maxi, mx, my, mz, mg = -1e90, -1, -1, -1, -1, -1
    if data[field].size > 0:
        maxi = na.argmax(data[field])
        ma = data[field][maxi]
        mx, my, mz = [data[ax][maxi] for ax in 'xyz']
        mg = data["GridIndices"][maxi]
    return (ma, maxi, mx, my, mz, mg)
def _combMaxLocation(data, *args):
    args = [na.atleast_1d(arg) for arg in args]
    i = na.argmax(args[0]) # ma is arg[0]
    return [arg[i] for arg in args]
add_quantity("MaxLocation", function=_MaxLocation,
             combine_function=_combMaxLocation, n_ret = 6)

def _MinLocation(data, field):
    """
    This function returns the location of the minimum of a set
    of fields.
    """
    ma, mini, mx, my, mz, mg = 1e90, -1, -1, -1, -1, -1
    if data[field].size > 0:
        mini = na.argmin(data[field])
        ma = data[field][mini]
        mx, my, mz = [data[ax][mini] for ax in 'xyz']
        mg = data["GridIndices"][mini]
    return (ma, mini, mx, my, mz, mg)
def _combMinLocation(data, *args):
    args = [na.atleast_1d(arg) for arg in args]
    i = na.argmin(args[0]) # ma is arg[0]
    return [arg[i] for arg in args]
add_quantity("MinLocation", function=_MinLocation,
             combine_function=_combMinLocation, n_ret = 6)


def _TotalQuantity(data, fields):
    """
    This function sums up a given field over the entire region

    :param fields: The fields to sum up
    """
    fields = ensure_list(fields)
    totals = []
    for field in fields:
        if data[field].size < 1:
            totals.append(0.0)
            continue
        totals.append(data[field].sum())
    return len(fields), totals
def _combTotalQuantity(data, n_fields, totals):
    totals = na.atleast_2d(totals)
    n_fields = totals.shape[1]
    return [na.sum(totals[:,i]) for i in range(n_fields)]
add_quantity("TotalQuantity", function=_TotalQuantity,
                combine_function=_combTotalQuantity, n_ret=2)
